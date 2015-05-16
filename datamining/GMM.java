import Jama.*;
import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

public class GMM {

    public static void main(String[] args) {
        GMM gmm = new GMM();
        gmm.test();
    }

    public void test() {
        GMMClustering gmm = new GMMClustering();
        ArrayList<ArrayList<Double>> data = readToPoints();
        ArrayList<ArrayList<Double>> re = gmm.cluster(data, data.size(), 3, 2);

        try {
            File file = new File("/Users/Alan/Documents/DM/dataset/clusters3-2.txt");
            if (!file.exists()) {
                file.createNewFile();
            }
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < re.size(); i++) {
                double max = Double.MIN_VALUE;
                int index = -1;
                for (int j = 0; j < re.get(i).size(); j++) {
                    if(re.get(i).get(j) > max){
                        max = re.get(i).get(j);
                        index = j;
                    }
                }

                if(index==0){
                    String w = data.get(i).get(0) + "\t" + data.get(i).get(1) + "\t" + "r";
                    bw.write(w);
                    bw.write("\n");
                }
                if(index==1){
                    String w = data.get(i).get(0) + "\t" + data.get(i).get(1) + "\t" + "g";
                    bw.write(w);
                    bw.write("\n");
                }
                if(index==2){
                    String w = data.get(i).get(0) + "\t" + data.get(i).get(1) + "\t" + "b";
                    bw.write(w);
                    bw.write("\n");
                }
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public ArrayList<ArrayList<Double>> readToPoints() {
        String path = "/Users/Alan/Documents/DM/dataset/dataset2.txt";
        ArrayList<ArrayList<Double>> points = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            String[] ws;
            while ((line = br.readLine()) != null) {
                ws = line.split("\t");
                ArrayList<Double> d = new ArrayList<>();
                d.add(Double.parseDouble(ws[0]));
                d.add(Double.parseDouble(ws[1]));
                points.add(d);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return points;
    }

    public class GMMClustering {

        private double oldExp = -9999999d;
        private final double iterationDifference = 0.000000000000000000000001;

        private ArrayList<ArrayList<Double>> Miu;
        private ArrayList<Double> Pai;
        private ArrayList<ArrayList<ArrayList<Double>>> Sig;

        public ArrayList<ArrayList<Double>> cluster(ArrayList<ArrayList<Double>> data, int num, int k, int dimen) {

            Miu = pickRandomPoints(data, num, k);
            Pai = new ArrayList<>();
            Sig = new ArrayList<>();

            int[] label = getLabel(data, Miu, num);

            int[] lebels = new int[k];
            for (int i = 0; i < num; i++) {
                lebels[label[i]]++;
            }

            for (int i = 0; i < k; i++) {
                Pai.add((double) (lebels[i]) / (double) (num));
            }

            for (int i = 0; i < k; i++) {
                ArrayList<ArrayList<Double>> oneLableData = new ArrayList<ArrayList<Double>>();
                for (int j = 0; j < num; j++) {
                    if (label[j] == i) {
                        oneLableData.add(data.get(j));
                    }
                }
                Sig.add(cov(oneLableData, dimen, num));
            }

            while (true) {
                ArrayList<ArrayList<Double>> pDensity = pdf(data, Miu, num, k, dimen);
                double[][] probablity = new double[num][k];
                for (int i = 0; i < num; i++) {
                    for (int j = 0; j < k; j++) {
                        probablity[i][j] = pDensity.get(i).get(j) * Pai.get(j);
                    }
                }

                double[] sumProb = sum(probablity, 2);
                for (int i = 0; i < num; i++) {
                    for (int j = 0; j < k; j++) {
                        probablity[i][j] = probablity[i][j] / sumProb[i];
                    }
                }

                double[] Nk = sum(probablity, 1);

                double[] NkBy1 = new double[Nk.length];
                for (int i = 0; i < Nk.length; i++) {
                    NkBy1[i] = 1 / Nk[i];
                }

                double[][] newMiuArr = new Matrix(diag(NkBy1)).times(new Matrix(probablity).transpose()).times(new Matrix(getArray(data))).getArray();

                Miu.clear();
                ArrayList<ArrayList<Double>> ll = new ArrayList<>(newMiuArr.length);
                for (double[] foo : newMiuArr) {
                    ArrayList<Double> l = new ArrayList<>();
                    for (double d : foo) {
                        l.add(d);
                    }
                    ll.add(l);
                }
                Miu = ll;

                double[][] newPieArr = new double[k][1];
                for (int i = 0; i < Nk.length; i++) {
                    newPieArr[i][0] = Nk[i] / num;
                }
                Pai.clear();
                Pai.add(newPieArr[0][0]);
                Pai.add(newPieArr[1][0]);

                Pai.add(newPieArr[2][0]);

                double[][][] newSig = new double[k][dimen][dimen];
                for (int i = 0; i < k; i++) {
                    double[][] xShift;
                    double[] rowmiu = newMiuArr[i];
                    double[][] mee = new double[num][2];
                    for (int mm = 0; mm < num; mm++) {
                        mee[mm] = rowmiu;
                    }

                    Matrix r = new Matrix(getArray(data));
                    Matrix rr = r.minus(new Matrix(mee));
                    xShift = rr.getArray();

                    double[] wRow = new double[num];
                    for (int j = 0; j < num; j++) {
                        wRow[j] = probablity[j][i];
                    }
                    double[][] diagWRow = diag(wRow);

                    double[][] SigK = new Matrix(xShift).transpose().times(new Matrix(diagWRow).times(new Matrix(xShift))).getArray();

                    for (int j = 0; j < dimen; j++) {
                        for (int l = 0; l < dimen; l++) {
                            newSig[i][j][l] = SigK[j][l] / Nk[i];
                        }
                    }
                }
                Sig.clear();
                for (int a = 0; a < newSig.length; a++) {
                    ArrayList<ArrayList<Double>> arrr = new ArrayList<>();
                    for (int b = 0; b < newSig[a].length; b++) {
                        double[] re = newSig[a][b];
                        List<Double> list = DoubleStream.of(re).boxed().collect(Collectors.toList());
                        ArrayList<Double> l = new ArrayList<Double>(list);
                        arrr.add(l);
                    }
                    Sig.add(arrr);
                }

                double[][] LikelyHood = new Matrix(getArray(pDensity)).times(new Matrix(newPieArr)).getArray();
                for (int i = 0; i < num; i++) {
                    LikelyHood[i][0] = Math.log(LikelyHood[i][0]);
                }
                double Exp = sum(LikelyHood, 1)[0];

                if (Exp - oldExp < iterationDifference) {
                    return pDensity;
                }
                oldExp = Exp;
            }
        }


        public ArrayList<ArrayList<Double>> pdf(ArrayList<ArrayList<Double>> data, ArrayList<ArrayList<Double>> Miu, int num, int k, int dimen) {
            double[][] px = new double[num][k];
            for (int i = 0; i < num; i++) {
                for (int j = 0; j < k; j++) {
                    ArrayList<Double> offset = minus(data.get(i), Miu.get(j));
                    ArrayList<ArrayList<Double>> s = Sig.get(j);

                    double[][] ss = getArray(Sig.get(j));

                    for (int p = 0; p < ss.length; p++) {
                        for (int q = 0; q < ss[0].length; q++) {
                           double v = ss[p][q];
                           ss[p][q] += Math.pow(1, -1);
                           s.get(p).set(q, v+ Math.pow(1, -32));
                           System.out.println(ss[p][q]);
                        }
                    }
                    double[][] inv = new Matrix(ss).inverse().getArray();

                    ArrayList<ArrayList<Double>> invSigma = new ArrayList<ArrayList<Double>>();
                    for (int m = 0; m < inv.length; m++) {
                        ArrayList<Double> arr = new ArrayList<Double>();
                        for (int n = 0; n < inv[m].length; n++) {
                            arr.add(inv[m][n]);
                        }
                        invSigma.add(arr);
                    }

                    double[][] ree = new Matrix(getOneDimenArray(offset))
                            .times(new Matrix(getArray(invSigma)))
                            .arrayTimes(new Matrix(getOneDimenArray(offset))).getArray();

                    double[] tmp = sum(ree, 2);
                    double coef = Math.pow((2 * Math.PI), -(double) dimen / 2d) * Math.sqrt(det(invSigma, invSigma.size()));
                    px[i][j] = coef * Math.pow(Math.E, -0.5 * tmp[0]);

                }
            }

            return getList(px);
        }

        public ArrayList<ArrayList<Double>> pickRandomPoints(ArrayList<ArrayList<Double>> data, int num, int k) {
            ArrayList<ArrayList<Double>> firstPoints = null;
            firstPoints = new ArrayList<>();
            Random rn = new Random();
            HashSet<Integer> set = new HashSet<>();
            for (int i = 0; i < k; i++) {
                int answer = rn.nextInt(num);
                if (!set.contains(answer)) {
                    set.add(answer);
                    firstPoints.add(data.get(answer));
                }
            }

            return firstPoints;
        }

        public int[] getLabel(ArrayList<ArrayList<Double>> data, ArrayList<ArrayList<Double>> Miu, int num) {
            int[] labels = new int[num];
            for (int i = 0; i < num; i++) {
                labels[i] = 0;

                double minDistance = distance(data.get(i), Miu.get(0));
                for (int j = 1; j < 3; j++) {
                    if (distance(data.get(i), Miu.get(j)) < minDistance) {
                        minDistance = distance(data.get(i), Miu.get(j));
                        labels[i] = j;
                    }
                }
            }
            return labels;
        }

        public ArrayList<ArrayList<Double>> cov(ArrayList<ArrayList<Double>> data, int dimen, int num) {
            ArrayList<ArrayList<Double>> res = new ArrayList<ArrayList<Double>>();

            double[] sum = new double[dimen];
            for (ArrayList<Double> dataRow : data) {
                for (int i = 0; i < dimen; i++) {
                    sum[i] += dataRow.get(i);
                }
            }
            for (int i = 0; i < dimen; i++) {
                sum[i] = sum[i] / num;
            }

            for (int i = 0; i < dimen; i++) {
                ArrayList<Double> tmp = new ArrayList<Double>();
                for (int j = 0; j < dimen; j++) {
                    double cov = 0;
                    for (ArrayList<Double> dataRow : data) {
                        cov += (dataRow.get(i) - sum[i]) * (dataRow.get(j) - sum[j]);
                    }
                    tmp.add(cov);
                }
                res.add(tmp);
            }
            return res;
        }

        public double[] sum(double[][] a, int mark) {
            double res[] = new double[a.length];
            if (mark == 1) {
                res = new double[a[0].length];
                for (int i = 0; i < a[0].length; i++) {
                    for (int j = 0; j < a.length; j++) {
                        res[i] += a[j][i];
                    }
                }
            } else if (mark == 2) {
                for (int i = 0; i < a.length; i++) {
                    for (int j = 0; j < a[0].length; j++) {
                        res[i] += a[i][j];
                    }
                }
            }
            return res;
        }

        public double[][] getArray(ArrayList<ArrayList<Double>> a) {
            int dimen = a.size();
            int dimen2 = a.get(0).size();

            double[][] res = new double[dimen][dimen2];

            for (int i = 0; i < dimen; i++) {
                for (int j = 0; j < dimen2; j++) {
                    res[i][j] = a.get(i).get(j);
                }
            }

            return res;
        }

        public double[][] diag(double[] a) {
            double[][] res = new double[a.length][a.length];
            for (int i = 0; i < a.length; i++) {
                for (int j = 0; j < a.length; j++) {
                    if (i == j) {
                        res[i][j] = a[i];
                    }
                }
            }
            return res;
        }

        public ArrayList<Double> minus(ArrayList<Double> a1, ArrayList<Double> a2) {
            ArrayList<Double> res = new ArrayList<Double>();
            for (int i = 0; i < a1.size(); i++) {
                res.add(a1.get(i) - a2.get(i));
            }
            return res;
        }

        public double[][] getOneDimenArray(ArrayList<Double> a) {
            int dimen = a.size();
            double[][] res = new double[1][dimen];

            for (int i = 0; i < dimen; i++) {
                res[0][i] = a.get(i);
            }

            return res;
        }

        public double det(ArrayList<ArrayList<Double>> data, int dimen) {
            return new Matrix(getArray(data)).det();
        }

        public ArrayList<ArrayList<Double>> getList(double[][] a) {
            ArrayList<ArrayList<Double>> res = new ArrayList<ArrayList<Double>>();
            for (int i = 0; i < a.length; i++) {
                ArrayList<Double> tmp = new ArrayList<Double>();
                for (int j = 0; j < a[i].length; j++) {
                    tmp.add(a[i][j]);
                }
                res.add(tmp);
            }
            return res;
        }

        public double distance(ArrayList<Double> d1, ArrayList<Double> d2) {
            double sum = 0;
            for (int i = 0; i < d1.size() - 1; i++) {
                sum += Math.pow(d1.get(i) - d2.get(i), 2);
            }
            return Math.sqrt(sum);
        }

    }
}
