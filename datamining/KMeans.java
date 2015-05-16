import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


public class KMeans {

    public static void main(String[] args) {
        KMeans km = new KMeans();
        km.cluster();
    }

    public void cluster() {

        double eps = 0.0000000000001;

        List<Point> points = readToPoints();

        int[] index = new int[3];
        Random randomGenerator = new Random();
        for (int idx = 0; idx < 3; idx++) {
            index[idx] = randomGenerator.nextInt(150);
        }

        Point centC1 = new Point(points.get(index[0]).getX(), points.get(index[0]).getY(), "centC1");
        Point centC2 = new Point(points.get(index[1]).getX(), points.get(index[1]).getY(), "centC2");
        Point centC3 = new Point(points.get(index[2]).getX(), points.get(index[2]).getY(), "centC3");

        List<Point> centers = new ArrayList<>(Arrays.asList(centC1, centC2, centC3));

        List<Point> cluster1 = new ArrayList<Point>();
        List<Point> cluster2 = new ArrayList<Point>();
        List<Point> cluster3 = new ArrayList<Point>();

        double diffC1 = Double.MAX_VALUE;
        double diffC2 = Double.MAX_VALUE;
        double diffC3 = Double.MAX_VALUE;

        while (diffC1 > eps || diffC2 > eps || diffC3 > eps) {
            cluster1.clear();
            cluster2.clear();
            cluster3.clear();

            for (Point p : points) {
                Point l = computeBelongTo(centers, p);
                if (l.getLabel() == "centC1") {
                    cluster1.add(p);
                }
                if (l.getLabel() == "centC2") {
                    cluster2.add(p);
                }
                if (l.getLabel() == "centC3") {
                    cluster3.add(p);
                }
            }

            Point newCentC1 = getCentroid(cluster1, "centC1");
            Point newCentC2 = getCentroid(cluster2, "centC2");
            Point newCentC3 = getCentroid(cluster3, "centC3");

            diffC1 = getDistance(newCentC1, centC1);
            diffC2 = getDistance(newCentC2, centC2);
            diffC3 = getDistance(newCentC3, centC3);

            centers.clear();
            centers.addAll(Arrays.asList(newCentC1, newCentC2, newCentC3));

            centC1 = newCentC1;
            centC2 = newCentC2;
            centC3 = newCentC3;

            System.out.println(diffC1 + " : " + diffC2 + " : " + diffC3);
        }

        List<List<Point>> clusters = new ArrayList<>(Arrays.asList(cluster1, cluster2, cluster3));

        try {
            File file = new File("/Users/Alan/Documents/DM/dataset/clusters1.txt");
            if (!file.exists()) {
                file.createNewFile();
            }
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            for (List<Point> ps : clusters) {
                for (Point p : ps) {
                    String w = p.getX() + "\t" + p.getY();
                    bw.write(w);
                    bw.write("\n");
                }
                bw.write("\n");
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Point computeBelongTo(List<Point> points, Point p) {
        double minDist = Double.MAX_VALUE;
        Point cp = null;
        for (Point pt : points) {
            if (getDistance(pt, p) < minDist) {
                minDist = getDistance(pt, p);
                cp = pt;
            }
        }
        return cp;
    }

    public List<Point> readToPoints() {
        String path = "/Users/Alan/Documents/DM/dataset/dataset1.txt";
        List<Point> points = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            String[] ws;
            while ((line = br.readLine()) != null) {
                ws = line.split("\t");
                Point p = new Point(Double.parseDouble(ws[0]), Double.parseDouble(ws[1]), ws[2]);
                points.add(p);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return points;
    }

    public Point getCentroid(List<Point> points, String l) {
        int num = points.size();
        double totalX = 0d;
        double totalY = 0d;

        for (int i = 0; i < points.size(); i++) {
            totalX += points.get(i).getX();
            totalY += points.get(i).getY();
        }
        return new Point(totalX / num, totalY / num, l);
    }

    public double getDistance(Point A, Point B) {

        return Math.sqrt(Math.pow((A.getX() - B.getX()), 2) + Math.pow((A.getY() - B.getY()), 2));
    }

    public class Point {
        private double X;
        private double Y;
        private String Label;

        public String getLabel() {
            return Label;
        }

        public void setLabel(String label) {
            Label = label;
        }

        public double getX() {
            return X;
        }

        public void setX(double x) {
            X = x;
        }

        public double getY() {
            return Y;
        }

        public void setY(double y) {
            Y = y;
        }

        public Point(double x, double y, String label) {
            X = x;
            Y = y;
            Label = label;
        }
    }
}
