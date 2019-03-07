/*
 * Data Hub Service (DHuS) - For Space data distribution.
 * Copyright (C) 2018 GAEL Systems
 *
 * This file is part of DHuS software sources.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
package fr.gael.drb.cortex.topic.sentinel1;

import static fr.gael.drb.cortex.topic.sentinel1.Sentinel1Utils.mercatorProjection;
import static fr.gael.drb.cortex.topic.sentinel1.Sentinel1Utils.mercatorToLatLon;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Objects;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;

/**
 * 2D concave hull algorithm.
 *
 * This is a Java implementation of "Algorithm 1" from the following article:
 *
 * "A new concave hull algorithm and concaveness measure for n-dimensional datasets."
 * Park, Jin-Seo, and Se-Jong Oh.
 * Journal of Information science and engineering 28.3 (2012): 587-600.
 */
public final class ConcaveHull2D
{
   private static final Logger LOGGER = Logger.getLogger(ConcaveHull2D.class);

   /**
    * Calls getConcaveHull, returns the convex hull if getConcaveHull throws an exception.
    *
    * @see #getConcaveHull(org.locationtech.jts.geom.GeometryFactory, org.locationtech.jts.geom.Coordinate[], double, double, double)
    *
    * @param geomFactory An instance of GeometryFactory (required)
    * @param points An array of points (required)
    * @param threshold See paper, 2.0 is often the best value
    * @param minLenght If the length of a segment is below that value, it will be skipped
    * @param minAngle Do not consider points that would make an angle larger than that value
    *
    * @return a concave hull or a convex hull if getConcaveHull failed
    */
   public static final Geometry getConcaveHullFailSafe(GeometryFactory geomFactory, Coordinate[] points,
         double threshold, double minLenght, double minAngle)
   {
      try
      {
         return getConcaveHull(geomFactory, points, threshold, minLenght, minAngle);
      }
      catch (Throwable t)
      {
         LOGGER.warn("Could not compute concave hull", t);
         return new ConvexHull(points, geomFactory).getConvexHull();
      }
   }

   /**
    * Computes the Concave hull of given array of points.
    *
    * @param geomFactory An instance of GeometryFactory (required)
    * @param points An array of points (required)
    * @param threshold See paper, 2.0 is often the best value
    * @param minLenght If the length of a segment is below that value, it will be skipped
    * @param minAngle Do not consider points that would make an angle larger than that value
    *
    * @return a concave hull
    *
    * @throws IllegalStateException If the length of the points parameter is lower than 3
    */
   @SuppressWarnings("unchecked")
   public static final Geometry getConcaveHull(GeometryFactory geomFactory, Coordinate[] points,
         double threshold, double minLenght, double minAngle)
   {
      // Check parameters
      Objects.requireNonNull(geomFactory, "GeometryFactory parameter must not be null");
      Objects.requireNonNull(points, "Coordinate[] parameter must not be null");
      if (points.length < 3)
      {
         throw new IllegalStateException("Not enough points (less than 3) to compute concave hull");
      }
      if (points.length == 3)
      {
         return geomFactory.createPolygon(points);
      }

      // List of points, ignore duplicate points
      CoordinateList pointList = new CoordinateList(points, false);

      // The convex hull is a subset of the concave hull that we want to compute
      Geometry convexHullGeom = new ConvexHull(pointList.toCoordinateArray(), geomFactory).getConvexHull();
      LinkedList<Coordinate> convexHull = new LinkedList<>(Arrays.asList(convexHullGeom.getCoordinates()));

      // Points to process (compute whether they are part of the concave hull or not)
      LinkedList<Coordinate> toProcess = new LinkedList<>(pointList);
      toProcess.removeAll(convexHull);

      // List of segments
      LinkedList<LineSegment> segments = makeSegmentList(convexHull.iterator());

      // Digging... for each segment
      ListIterator<LineSegment> it = segments.listIterator();
      LineSegment previous;
      LineSegment current = segments.peekLast();
      LineSegment next = it.next();
      while (it.hasNext()) // For each segment in order
      {
         previous = current;
         current  = next;
         next     = it.next();

         if (LOGGER.isDebugEnabled())
         {
            LOGGER.debug("Segment to test: " + current.toGeometry(geomFactory).toText());
         }

         double eh = current.getLength();
         if (eh < minLenght) continue; // Do not dig further if segment is below given threshold length

         Coordinate point = findBestCandidate(current, previous, next, toProcess, minAngle);

         if (point == null) continue; // No point found

         double dd = Math.min(point.distance(current.p0), point.distance(current.p1));
         if (eh/dd <= threshold) continue; // If point has not correct concaveness

         LineSegment segA = new LineSegment(current.p0, point);
         LineSegment segB = new LineSegment(point, current.p1);

         if (LOGGER.isDebugEnabled())
         {
            LOGGER.debug(current.toGeometry(geomFactory).toText()
                  + " -> " + segA.toGeometry(geomFactory).toText() + ", " + segB.toGeometry(geomFactory).toText());
         }

         // Remove current segment [Sa, Sb]
         it.previous();it.previous(); // previous must be called twice!!
         it.remove();

         it.add(segA); // add new segment Sa->P
         it.add(segB); // add new segment P->Sb
         it.previous(); // rewind so next() is segB

         // Set current and next for the next iteration
         current = previous;
         next = segA;

         // Remove processed point
         toProcess.remove(point);

         if (LOGGER.isDebugEnabled())
         {
            LOGGER.debug("New hull: " + geomFactory.createPolygon(toCoordinateList(segments).toCoordinateArray()).toText());
         }
      }

      // Creates result polygon
      CoordinateList result = toCoordinateList(segments);
      try
      {
         return geomFactory.createPolygon(result.toCoordinateArray());
      }
      catch (IllegalArgumentException ex)
      {
         LOGGER.warn("Could not create polygon from points: "
               + geomFactory.createMultiPointFromCoords(result.toCoordinateArray()).toText());
         throw ex;
      }
   }

   private static CoordinateList toCoordinateList(LinkedList<LineSegment> segments)
   {
      CoordinateList result = new CoordinateList();
      segments.stream()
              .<Coordinate>map((t) -> t.p0)
              .forEach((t) -> result.add(t, false));
      result.add(segments.peek().p0);
      return result;
   }

   private static Coordinate findBestCandidate(LineSegment segment, LineSegment previous, LineSegment next, List<Coordinate> points, double minAngle)
   {
      Coordinate res = null;
      double resDist = Double.MAX_VALUE;

      for (Coordinate p: points)
      {
         double distToSegment = pointToSegment(p, segment.p0, segment.p1);
         double distToPrevious = pointToSegment(p, previous.p0, previous.p1);
         double distToNext = pointToSegment(p, next.p0, next.p1);

         if (distToSegment < resDist && distToSegment < distToNext && distToSegment < distToPrevious)
         {
            double angle = Angle.angleBetween(segment.p0, p, segment.p1);
            if (angle < minAngle) continue; // angle threshold

            res = p;
            resDist = distToSegment;
         }
      }

      return res;
   }

   private static LinkedList<LineSegment> makeSegmentList(Iterator<Coordinate> input)
   {
      LinkedList<LineSegment> res = new LinkedList<>();
      Coordinate first = input.next();
      Coordinate previous = first;
      while (input.hasNext())
      {
         Coordinate current = input.next();
         res.add(new LineSegment(previous, current));
         previous = current;
      }
      res.add(new LineSegment(previous, first));
      return res;
   }

   // More robust than Distance.pointToSegment(...)
   private static double pointToSegment(Coordinate p, Coordinate A, Coordinate B)
   {
      // if start = end, then just compute distance to one of the endpoints
      if (A.x == B.x && A.y == B.y)
      {
         return p.distance(A);
      }

      // r_pXY = Xp dot XY
      double r_pAB = ((p.x - A.x) * (B.x - A.x) + (p.y - A.y) * (B.y - A.y));
      double r_pBA = ((p.x - B.x) * (A.x - B.x) + (p.y - B.y) * (A.y - B.y));

      if (r_pAB <= 0.0)
      {
         return p.distance(A);
      }
      if (r_pBA <= 0.0)
      {
         return p.distance(B);
      }

      // Distance between P and its orthogonal projection on [AB]
      double len2 = (B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y);
      double s = ((A.y - p.y) * (B.x - A.x) - (A.x - p.x) * (B.y - A.y)) / len2;
      return Math.abs(s) * Math.sqrt(len2);
   }

   /**
    * To test the concave hull algorithm.
    *
    * @param args args[0] must be a WKT encoded MULTIPOINT
    *
    * @throws ParseException malformed WKT string
    */
   public static void main(String[] args) throws ParseException
   {
      LOGGER.setLevel(Level.DEBUG);
      GeometryFactory geomFactory = new GeometryFactory();

      Geometry multipoint = new WKTReader().read(args[0]);
      Coordinate[] points = multipoint.getCoordinates();

      Coordinate[] projPoints = mercatorProjection(points);

      Geometry concaveHull = getConcaveHull(geomFactory, projPoints, 2., 3., Math.PI/2.);

      mercatorToLatLon(points, concaveHull.getCoordinates());
      concaveHull.geometryChanged();

      System.out.println(concaveHull.toText());
   }

}
