package neqsim.thermo.util.referenceequations;

/**
 * Simplified IAPWS-IF97 implementation providing a few basic water and steam
 * properties. Units follow the IAPWS convention: pressure in MPa, temperature
 * in Kelvin. Enthalpy is returned in kJ/kg and specific volume in m^3/kg.
 *
 * <p>This class only covers region 1 and 2 of the IAPWS industrial
 * formulation and is not intended for highly accurate work, but gives
 * reasonable values for common engineering calculations.</p>
 */
public final class IAPWSIF97 {
  private IAPWSIF97() {}

  private static final double[] I1 = new double[] {0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32};
  private static final double[] J1 = new double[] {-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41};
  private static final double[] N1 = new double[] {
    0.14632971213167,-0.84548187169114,-3.7563603672040,3.3855169168385,
    -0.95791963387872,0.15772038513228,-0.016616417199501,0.00081214629983568,
    0.00028319080123804,-0.00060706301565874,-0.018990068218419,
    -0.032529748770505,-0.021841717175414,-5.283835796993e-05,
    -0.00047184321073267,-0.00030001780793026,4.7661393906987e-05,
    -4.4141845330846e-06,-7.2694996297594e-16,-3.1679644845054e-05,
    -2.8270797985312e-06,-8.5205128120103e-10,-2.2425281908e-06,
    -6.5171222895601e-07,-1.4341729937924e-13,-4.0516996860117e-07,
    -1.2734301741641e-09,-1.7424871230634e-10,-6.8762131295531e-19,
    1.4478307828521e-20,2.6335781662795e-23,-1.1947622640071e-23,
    1.8228094581404e-24,-9.3537087292458e-26};

  private static final double[] J0 = new double[] {0,1,-5,-4,-3,-2,-1,2,3};
  private static final double[] N0 = new double[] {
    -9.6927686500217,10.086655968018, -0.005608791128302,0.071452738081455,
    -0.40710498223928,1.4240819171444,-4.383951131945,-0.28408632460772,
    0.021268463753307};

  private static final double[] IR2 = new double[] {1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,8,8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24};
  private static final double[] JR2 = new double[] {0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,25,8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58};
  private static final double[] NR2 = new double[] {
    -0.0017731742473213,-0.017834862292358,-0.045996013696365,-0.057581259083432,
    -0.05032527872793,-3.3032641670203e-05,-0.00018948987516315,
    -0.0039392777243355,-0.043797295650573,-2.6674547914087e-05,
    2.0481737692309e-08,4.3870667284435e-07,-3.227767723857e-05,
    -0.0015033924542148,-0.040668253562649,-7.8847309559367e-10,
    1.2790717852285e-08,4.8225372718507e-07,2.2922076337661e-06,
    -1.6714766451061e-11,-0.0021171472321355,-23.895741934104,
    -5.905956432427e-18,-1.262180889911e-06,-0.038946842435739,
    1.1256211360459e-11,-8.2311340897998,1.9809712802088e-08,
    1.0406965210174e-19,-1.0234747095929e-13,-1.0018179379511e-09,
    -8.0882908646985e-11,0.10693031879409,-0.33662250574171,
    8.9185845355421e-25,3.0629316876232e-13,-4.2002467698208e-06,
    -5.9056029685639e-26,3.7826947613457e-06,-1.2768608934681e-15,
    7.3087610595061e-29,5.5414715350778e-17,-9.436970724121e-07};

  private static final double R = 0.461526; // kJ/(kg K)

  /** Enthalpy h(p,T) in kJ/kg. */
  public static double h_pT(double pMPa, double temperatureK) {
    if (temperatureK <= 623.15 && pMPa >= psat_T(temperatureK)) {
      return R * 1386.0 * dGammaTau1(pMPa, temperatureK);
    }
    return R * 540.0 * dGammaTau2(pMPa, temperatureK);
  }

  /** Specific volume v(p,T) in m^3/kg. */
  public static double v_pT(double pMPa, double temperatureK) {
    if (temperatureK <= 623.15 && pMPa >= psat_T(temperatureK)) {
      return 1e-3 * R * temperatureK / 16.53 * dGammaPi1(pMPa, temperatureK);
    }
    return 1e-3 * R * temperatureK * dGammaPi2(pMPa, temperatureK);
  }

  /** Saturation pressure as function of temperature. */
  public static double psat_T(double temperatureK) {
    double Tmin = 273.16;
    double Tc = 647.096;
    if (temperatureK < Tmin || temperatureK > Tc) {
      return Double.NaN;
    }
    double[] n = {1167.05214527670,-724213.167032060,-17.0738469400920,
        12020.8247024700,-3232555.03223330,14.9151086135300,
        -4823.26573615910,405113.405420570,-0.238555575678490,
        650.175348447980};
    double upsilon = temperatureK + n[8] / (temperatureK - n[9]);
    double A = (upsilon + n[0]) * upsilon + n[1];
    double B = (n[2] * upsilon + n[3]) * upsilon + n[4];
    double C = (n[5] * upsilon + n[6]) * upsilon + n[7];
    double beta = 2 * C / (-B + Math.sqrt(B * B - 4 * A * C));
    return Math.pow(beta, 4);
  }

  /** Saturation temperature as function of pressure. */
  public static double Tsat_p(double pressureMPa) {
    double Tmin = 273.16;
    double pmin = psat_T(Tmin);
    double pc = 22.064;
    if (pressureMPa < pmin || pressureMPa > pc) {
      return Double.NaN;
    }
    double[] n = {1167.05214527670,-724213.167032060,-17.0738469400920,
        12020.8247024700,-3232555.03223330,14.9151086135300,
        -4823.26573615910,405113.405420570,-0.238555575678490,
        650.175348447980};
    double beta = Math.pow(pressureMPa, 0.25);
    double E = (beta + n[2]) * beta + n[5];
    double F = (n[0] * beta + n[3]) * beta + n[6];
    double G = (n[1] * beta + n[4]) * beta + n[7];
    double D = 2 * G / (-F - Math.sqrt(F * F - 4 * E * G));
    double theta = (n[9] + D - Math.sqrt(Math.pow(n[9] + D, 2)
        - 4 * (n[8] + n[9] * D))) / 2.0;
    return theta;
  }

  private static double dGammaTau1(double p, double T) {
    double pi = p / 16.53;
    double tau = 1386.0 / T;
    double sum = 0.0;
    for (int i = 0; i < N1.length; i++) {
      sum += N1[i] * Math.pow(7.1 - pi, I1[i]) * J1[i]
          * Math.pow(tau - 1.222, J1[i] - 1.0);
    }
    return sum;
  }

  private static double dGammaPi1(double p, double T) {
    double pi = p / 16.53;
    double tau = 1386.0 / T;
    double sum = 0.0;
    for (int i = 0; i < N1.length; i++) {
      sum += -N1[i] * I1[i] * Math.pow(7.1 - pi, I1[i] - 1.0)
          * Math.pow(tau - 1.222, J1[i]);
    }
    return sum;
  }

  private static double dGammaTau2(double p, double T) {
    double pi = p;
    double tau = 540.0 / T;
    double sum0 = 0.0;
    for (int i = 0; i < N0.length; i++) {
      sum0 += N0[i] * J0[i] * Math.pow(tau, J0[i] - 1.0);
    }
    double sumR = 0.0;
    for (int i = 0; i < NR2.length; i++) {
      sumR += NR2[i] * Math.pow(pi, IR2[i]) * JR2[i]
          * Math.pow(tau - 0.5, JR2[i] - 1.0);
    }
    return sum0 + sumR;
  }

  private static double dGammaPi2(double p, double T) {
    double pi = p;
    double tau = 540.0 / T;
    double sum0 = 1.0 / pi;
    double sumR = 0.0;
    for (int i = 0; i < NR2.length; i++) {
      sumR += NR2[i] * IR2[i] * Math.pow(pi, IR2[i] - 1.0)
          * Math.pow(tau - 0.5, JR2[i]);
    }
    return sum0 + sumR;
  }
}
