/*
 * Data Hub Service (DHuS) - For Space data distribution.
 * Copyright (C) 2014,2015,2018 GAEL Systems
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

import fr.gael.drb.DrbNode;
import fr.gael.drb.DrbSequence;
import fr.gael.drb.query.Query;
import fr.gael.drbx.image.DrbCollectionImage;

import java.util.Scanner;

import javax.media.jai.RenderedImageAdapter;

import org.apache.log4j.Logger;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;

public class Sentinel1Utils
{
   static private final Logger LOGGER = Logger.getLogger(Sentinel1Utils.class);

   private Sentinel1Utils() {}

   public static boolean isAscending(Object source)
   {
      Object wrapped = source;

      if (source instanceof RenderedImageAdapter)
      {
         wrapped = ((RenderedImageAdapter) source).getWrappedImage();
      }

      if (wrapped instanceof DrbCollectionImage)
      {
         DrbNode node = (DrbNode) ((DrbCollectionImage) wrapped).getItemSource();

         Query q = new Query(
               "data("
               + "manifest.safe/XFDU/metadataSection/"
               + "metadataObject[fn:matches(@ID,\".+OrbitReference\")]/"
               + "metadataWrap/xmlData/"
               + "orbitReference/extension/orbitProperties/pass)");

         DrbSequence seq = q.evaluate(node);

         if ((seq == null) || (seq.getLength() < 1))
         {
            LOGGER.error("Cannot extract ASC/DESC pass from the product");
            return false;
         }

         String mode = seq.getItem(0).getValue().toString();
         if (mode.toLowerCase().equals("ascending"))
         {
            return true;
         }
         else if (mode.toLowerCase().equals("descending"))
         {
            return false;
         }
         else
         {
            LOGGER.error("Unknown mode : " + mode);
         }
      }
      LOGGER.error("Pass mode connot be extracted");
      return false;
   }

   public static Coordinate[] mercatorProjection(Coordinate[] points)
   {
      double toRad = Math.PI / 180.;
      double toDeg = 180. / Math.PI;

      Coordinate[] res = new Coordinate[points.length];
      for (int it=0; it<res.length; it++)
      {
         double lat = points[it].y * toRad;
         res[it] = new Coordinate(
               points[it].x,
               Math.log( Math.tan(lat) + 1./Math.cos(lat) ) * toDeg,
               it); // HACK: uses the Z coordinate to store the index in the source array to recover the original value
      }
      return res;
   }

   public static void mercatorToLatLon(Coordinate[] source, Coordinate[] data)
   {
      for (Coordinate point: data)
      {
         int index = (int)point.z; // Uses the HACK from #mercatorProjection(Coordinate[])
         point.x = source[index].x;
         point.y = source[index].y;
      }
   }

   public static String makeConcaveHull(String points)
   {
      Scanner scanner = new Scanner(points);
      scanner.useDelimiter(" ");

      CoordinateList coordList = new CoordinateList();

      // Parse point, check for antemeridian crossing
      boolean fixTopo = false;
      Coordinate prev = null;
      while (scanner.hasNext())
      {
         String[] components = scanner.next().split(",");
         // Lat,Lon --> Lon,Lat (x,y)
         double x = Double.parseDouble(components[1]);
         double y = Double.parseDouble(components[0]);

         Coordinate elem = new Coordinate(x, y);
         coordList.add(elem, false);

         if (!fixTopo && prev != null && Math.abs(x - prev.x) > 330.)
         {
            LOGGER.debug("Two adjacent points seem to be crossing the ante-meridian line, topology has to be fixed");
            fixTopo = true;
         }
         prev = elem;
      }
      if (fixTopo)
      {
         // Fix topology
         for (Object coordObj: coordList)
         {
            Coordinate coord = (Coordinate)coordObj;
            if (coord.x < 0) coord.x = coord.x + 360.0;
         }
      }

      GeometryFactory geomFactory = new GeometryFactory();

      // Debugging
      if (LOGGER.isDebugEnabled())
      {
         LOGGER.debug("makeConcaveHull input = " + geomFactory.createMultiPointFromCoords(coordList.toCoordinateArray()).toText());
      }

      // Compute concave hull
      Coordinate[] source = coordList.toCoordinateArray();
      Coordinate[] projPoints = mercatorProjection(source);
      Geometry enveloppe = ConcaveHull2D.getConcaveHullFailSafe(geomFactory, projPoints, 2., 3., Math.PI/2.);
      mercatorToLatLon(source, enveloppe.getCoordinates());
      enveloppe.geometryChanged();

      // Debugging
      if (LOGGER.isDebugEnabled())
      {
         LOGGER.debug("makeConcaveHull output = " + enveloppe.toText());
      }

      // output
      StringBuilder sb = new StringBuilder(enveloppe.getNumPoints() * 2 * 5);

      Coordinate[] res = enveloppe.getCoordinates();
      for (int it=0; it<res.length; it++)
      {
         double lat = res[it].y;
         double lon = res[it].x;

         // Fix topology
         if (fixTopo && lon > 180) lon = lon - 360.0;

         sb.append(Double.toString(lat)).append(',').append(Double.toString(lon)).append(' ');
      }

      return sb.toString();
   }
}
