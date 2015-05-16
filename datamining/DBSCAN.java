import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

// implementation of DBSCAN algorithm 
public class DBSCAN {

    public static void main(String[] args) {
        DBSCAN dbscan = new DBSCAN();
        dbscan.DBScan();
    }

    public void DBScan() {
        double eps = 0.9d;
        int minpts = 5;
        List<Point> points = readToPoints();
        List<List<Point>> allClusters = new ArrayList<>();

        for (Point p : points) {
            if (p.isVisited()) {
                continue;
            }
            List<Point> clusterPoints = new ArrayList<Point>();

            points.stream().filter(pt -> getDistance(pt, p) <= eps).forEach(clusterPoints::add);

            if (clusterPoints.size() >= minpts) {
                p.setCore(true);
                p.setVisited(true);
                if (!clusterPoints.isEmpty()) {
                    allClusters.add(clusterPoints);
                }
            }
        }

        for (int i = 0; i < allClusters.size(); i++) {
            for (int j = 0; j < allClusters.size(); j++) {
                if (i != j) {
                    if (allClusters.get(i) == null || allClusters.get(j) == null) {
                        continue;
                    }
                    for (int k = 0; k < allClusters.get(j).size(); k++) {
                        Point p = allClusters.get(j).get(k);
                        if (p.isCore() && allClusters.get(i).contains(p)) {
                            for (int l = 0; l < allClusters.get(j).size(); l++) {
                                if (!allClusters.get(i).contains(allClusters.get(j).get(l))) {
                                    allClusters.get(i).add(allClusters.get(j).get(l));
                                }
                            }
                            allClusters.get(j).clear();
                            break;
                        }
                    }
                }
            }
        }

        List<List<Point>> result = allClusters.stream().filter(ls -> !ls.isEmpty()).collect(Collectors.toCollection(ArrayList::new));

        try {
            File file = new File("/Users/Alan/Documents/DM/dataset/clusters2.txt");
            if (!file.exists()) {
                file.createNewFile();
            }
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            for (List<Point> ps : result) {
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

    public class Point {
        private double X;
        private double Y;
        private boolean isCore;
        private boolean isVisited;

        public Point(double x, double y, boolean iscore, boolean isvisited) {
            X = x;
            Y = y;
            isCore = iscore;
            isVisited = isvisited;
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

        public boolean isCore() {
            return isCore;
        }

        public void setCore(boolean isCore) {
            this.isCore = isCore;
        }

        public boolean isVisited() {
            return isVisited;
        }

        public void setVisited(boolean isVisited) {
            this.isVisited = isVisited;
        }
    }

    public double getDistance(Point A, Point B) {
        return Math.sqrt(Math.pow((A.getX() - B.getX()), 2) + Math.pow((A.getY() - B.getY()), 2));
    }

    public List<Point> readToPoints() {
        String path = "/Users/Alan/Documents/DM/dataset/dataset2.txt";
        List<Point> points = new ArrayList<Point>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            String[] ws;
            while ((line = br.readLine()) != null) {
                ws = line.split("\t");
                Point p = new Point(Double.parseDouble(ws[0]), Double.parseDouble(ws[1]), false, false);
                points.add(p);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return points;
    }

}
