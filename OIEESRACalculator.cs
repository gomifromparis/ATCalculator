using System;

namespace OIECalculators
{
    //ORIGINAL C CODE IS FROM OIE, Ecole des Mines de Paris, see http://www.oie.mines-paristech.fr/Valorisation/Outils/Clear-Sky-Library/
    // I merely translated it to C#

    public static class ESRACalculator
    {

        /*-------------------------------------------------------------------------*/
        /*                        ECOLE DES MINES DE PARIS                         */
        /*        CENTRE D'ENERGETIQUE - GROUPE TELEDETECTION & MODELISATION       */
        /*                       Rue Claude Daunesse, BP 207                       */
        /*                   06904 Sophia Antipolis cedex, FRANCE                  */
        /*          Tel (+33) 04 93 95 74 49     Fax (+33) 04 93 95 75 35          */
        /*                       E-mail : (name)@cenerg.cma.fr                     */
        /*-------------------------------------------------------------------------*/
        /* csmodels_lib.c                                                          */
        /* O. Bauer - February 1997 modified M. Lefevre December 2001              */
        /* modified L. Wald September 2014 (bug correction in irradiance functions)*/
        /* Some useful functions for clear sky models for solar irradiation and    */
        /* irradiance computations                                                 */
        /* cf. Rigollier C., Bauer O., Wald L., 2000. On the clear sky model of the*/
        /* 4th european solar radiation atlas with respect to the Heliosat method  */
        /* Solar Energy, 68(1), 33-48p.                                            */
        /*-------------------------------------------------------------------------*/


        /* NOTATIONS */
        /*************/
        /* phi    : latitude of the site, positive to North */
        /* lambda : longitude of the site, positive to East */
        /* delta  : solar declination angle */
        /* omega  : solar hour angle */
        /* gamma  : solar altitude (or elevation) angle */
        /* theta  : solar incidence (or zenithal) angle */
        /* alpha  : solar azimuthal angle */

        /* t   : solar time = true solar time (TST) = local apparent time (LAT) */
        /* LAT : local apparent time or solar time or true solar time (TST)
           --> this system of time offers the advantage of symmetry of the solar 
               geometry about the north-south line */
        /* LMT : local mean time or clock time */
        /* UT  : Universal Time, is GMT measured from Greenwich mean midnight */

        /* omega_sr : sunrise hour angle */
        /* omega_ss : sunset hour angle */
        /* t_sr     : time of astronomical sunrise */
        /* t_ss     : time of astronomical sunset */

        /* omega1 : solar hour angle at beginning of the time period */
        /* omega2 : solar hour angle at end of the time period */

        /* S0  : astronomical daylength or astronomical sunshine duration */
        /* I0  : solar constant = annual mean value of extraterrestrial direct solar 
           irradiance G0 (1367.0 W/m2) */
        /* G0  : extraterrestrial global solar radiation (on an horizontal plane) =B0*/
        /* G0h : hourly extraterrestrial solar irradiation (on an horizontal plane) */
        /* G0d : daily extraterrestrial solar irradiation (on an horizontal plane) */

        /* z      : station height above sea level */
        /* p,p0   : pressure at z, and at sea level */
        /* height_corr : p/p0 */
        /* m      : relative optical air mass */
        /* deltaR : integral Rayleigh optical thickness */
        /* TL     : Linke turbidity factor */
        /* TLAM2  : Linke turbidity factor for air mass 2 i.e. TLAM2 = TL(2) */

        /* Gc : cloudless (or clear) sky global horizontal irradiance/irradiation 
              = global irradiance under clear sky conditions on an horizontal surface*/
        /* Bc': cloudless (or clear) sky direct beam irradiance */
        /* Bc : cloudless (or clear) sky beam horizontal irradiance/irradiation 
           Bc = Bc'.Math.Sin(gamma) */
        /* Dc : cloudless (or clear) sky diffuse horizontal irradiance/irradiation */
        /* Gch, Bch, Dch : integrated hourly values of global, beam and diffuse 
           horizontal irradiation */
        /* Gcd, Bcd, Dcd : integrated daily values of global, beam and diffuse 
           horizontal irradiation */
        /* Trb : transmittance at zenith for beam radiation */
        /* Trd : transmittance at zenith for diffuse radiation */
        /* Fb  : beam angular function */
        /* Fd  : diffuse angular function */

        /* NB : All angles are computed in radians as a standard.
                The basic trigonometric calculations on the position of the sun are 
            carried out in LAT. */


        /***********************************************************************/
        /*                                                                     */
        /* O P T I C A L  P R O P E R T I E S  O F  T H E  A T M O S P H E R E */
        /*                                                                     */
        /***********************************************************************/

        /***********************************/
        /* STATION HEIGHT CORRECTION, p/p0 */
        /***********************************/

        public static int station_height_correction(double z, out double height_corr)
        /* Source : Rodgers, Souster and Page (1978) */
        /* Inputs :
           z : station height above sea level (in meters) */
        /* Outputs :
           height_corr : station height correction */
        /* The procedure "station_height_correction" computes the station height 
           correction i.e. the ratio of mean atmospheric pressure (p) at the site 
           elevation, to the mean atmospheric pressure at sea level (p0). This ratio is
           used to correct the optical path length to take account of station height, 
           and is important in mountainous areas. Returns 0 if OK, 1 otherwise. */
        {
            if (z < 4000.0)
                height_corr = 1.0 - (z / 1.0e+04);
            else
                height_corr = Math.Exp(-(0.1174 + (0.0017 * z / 1.0e+03)) * z / 1.0e+03);
            return 0;
        }

        public static int station_height_correction2(double z, out double height_corr)
        /* Source : */
        /* Inputs :
           z : station height above sea level (in meters) */
        /* Outputs :
           height_corr : station height correction */
        /* The procedure "station_height_correction" computes the station height 
           correction i.e. the ratio of mean atmospheric pressure (p) at the site 
           elevation, to the mean atmospheric pressure at sea level (p0). This ratio is
           used to correct the optical path length to take account of station height, 
           and is important in mountainous areas. Returns 0 if OK, 1 otherwise. */
        {
            const double zh = 8434.5; /* scale height of the Rayleigh atmosphere (in 
			      meters) in the Standard Atmosphere */
            height_corr = Math.Exp(-z / zh);
            return 0;
        }

        /********************************/
        /* RELATIVE OPTICAL AIR MASS, m */
        /********************************/


        public static int air_mass_Kasten(double gamma, double z, out double m)

        /* Source :
           Kasten, F. and Young, A.T., 1989.
           Revised optical air mass tables and approximation formula. Applied Optics, 
           28 (22), 4735-4738 */
        /* Inputs :
           gamma  : solar elevation angle (in radians)
           z      : station height above sea level (in meters) */
        /* Outputs :
           m : relative optical air mass */
        /* The procedure "air_mass_Kasten" provides the relative optical air mass.
           It expresses the optical path length of the solar beam through the 
           atmosphere as a ratio to sea level optical path through a standard 
           atmosphere with the sun at the zenith. Returns 0 if OK, 1 otherwise. */
        /* Remarks:
           - if gamma=90 deg then gamma_cor=0.010031 deg, gamma(true)=90.010031 deg 
           and m = 1.000 ;
           - if gamma=0 deg then gamma_corr=0.56038851 deg, gamma(true)=0.56038851 deg
           and m= 30.666; without taking account of the refraction m= 37.920 */
        {
            const double deg_rad = (Math.PI / 180.0); /* converts decimal degrees into radians*/
            const double rad_deg = (180.0 / Math.PI); /* converts radians into decimal degrees*/
            double height_corr, gamma_corr, gamma_deg;
            double c1, c2, c3, c4, c5, c6, a, b, c;
            int ier;

            ier = 1;
            m = Double.NaN;
            /* Correction of solar elevation for the effect of atmospheric refraction */
            c1 = 0.061359;
            c2 = 0.1594;
            c3 = 1.1230;
            c4 = 0.065656;
            c5 = 28.9344;
            c6 = 277.3971;
            /* if(gamma <= 0.0) */  /* Aguiar */
            if (gamma < 0.0)
                m = 0.0;
            else
            {
                gamma_corr = c1 * (c2 + c3 * gamma + c4 * gamma * gamma) /
                                  (1.0 + c5 * gamma + c6 * gamma * gamma);
                gamma = gamma + gamma_corr;
            }

            gamma_deg = gamma * rad_deg;
            a = 0.50572;
            b = 6.07995;  /* in degrees */
            c = -1.6364;
            ier = station_height_correction2(z, out height_corr);
            if ((ier == 0) && (gamma > 0.0))
            {
                ier = 0;
                m = height_corr / (Math.Sin(gamma) + a * Math.Pow((gamma_deg + b), c));
            }

            return (ier);
        }

        /*************************************************************/
        /* INTEGRAL RAYLEIGH OPTICAL THICKNESS PER UNIT MASS, deltaR */
        /*************************************************************/

        public static int deltaR_Kasten_old(double m, out double deltaR)

        /* Source :
           Kasten, F., 1980.
           A simple parameterization of the pyrheliometric formula for determining the
           Linke turbidity factor. Meteorol. Rundschau, 33, 124-127.*/
        /* Inputs :
           m : relative optical air mass */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Kasten_old" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            deltaR = Double.NaN;
            if (m <= 0.0)
                ier = 1;
            else
            {
                ier = 0;
                deltaR = 1.0 / (9.4 + 0.9 * m);
            }
            return (ier);
        }

        public static int deltaR_Louche(double m, out double deltaR)

        /* Source :
           Louche, A., Peri, G. and Iqbal, M., 1986.
           An analysis of Linke turbidity factor. Solar Energy, 37 (6), 393-396. */
        /* Inputs :
           m : relative optical air mass */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Louche" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            double a, b, c, d, e, m2, m3, m4;

            ier = 0;
            a = 6.5567;
            b = 1.7513;
            c = -0.1202;
            d = 0.0065;
            e = -0.00013;
            /* if( (m >= 1.0) && (m <= 20.0) ) */
            m2 = m * m;
            m3 = m2 * m;
            m4 = m3 * m;
            deltaR = 1.0 / (a + b * m + c * m2 + d * m3 + e * m4);

            return (ier);
        }

        public static int deltaR_Grenier(double m, out double deltaR)

        /* Source :
           Grenier, J.C., De La Casiniere, A. and Cabot, T., 1994.
           A spectral model of Linke's turbidity factor and its experimental 
           implications. Solar Energy, 52 (4), 303-313. */
        /* Inputs :
           m : relative optical air mass */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Grenier" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            double a, b, c, d, e, m2, m3, m4;

            ier = 0;
            a = 5.4729;
            b = 3.0312;
            c = -0.6329;
            d = 0.0910;
            e = -0.00512;
            /* if( (m >= 1.0) && (m <= 6.0) ) */
            m2 = m * m;
            m3 = m2 * m;
            m4 = m3 * m;
            deltaR = 1.0 / (a + b * m + c * m2 + d * m3 + e * m4);

            return (ier);
        }

        public static int deltaR_Kasten(double m, out double deltaR)

        /* Source :
           Kasten, F., 1996. (adjusted Louche et al., 1986)
           The Linke turbidity factor based on improved values of the integral Rayleigh
           optical thickness. Solar Energy, 56 (3), 239-244. */
        /* Inputs :
           m : relative optical air mass */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Kasten" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            double a, b, c, d, e, m2, m3, m4;

            ier = 0;
            a = 6.6296;
            b = 1.7513;
            c = -0.1202;
            d = 0.0065;
            e = -0.00013;
            /* if(Math.Abs(m-1.0) < 0.0) */
            if (m <= 0.0)
                deltaR = 0.0;
            else
            {
                m2 = m * m;
                m3 = m2 * m;
                m4 = m3 * m;
                deltaR = 1.0 / (a + b * m + c * m2 + d * m3 + e * m4);
            }

            return (ier);
        }

        public static int deltaR_Page_old(double m, out double deltaR)

        /* Source :
           cf. Rigollier et al. 2000 */
        /* Inputs :
           m : relative optical air mass */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Page" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            double a, b, c, d, e, m2, m3, m4;

            ier = 0;
            a = 6.6296;
            b = 1.7513;
            c = -0.1202;
            d = 0.0065;
            e = -0.00013;
            if (m <= 0.0)
                deltaR = 0.0;
            else
            {
                if (m <= 20)
                {
                    m2 = m * m;
                    m3 = m2 * m;
                    m4 = m3 * m;
                    deltaR = 1.0 / (a + b * m + c * m2 + d * m3 + e * m4);
                }
                else  /* modification by Page (1996) */
                    deltaR = 1.0 / (10.4 + 0.718 * m);
            }

            return (ier);
        }

        public static int deltaR_Page(double msl, double m, double height_corr, out double deltaR)

        /* Source :
           Page, J., 2001. (modified Rigollier et al. 2000 for pressure correction) */
        /* Inputs :
           m   : relative optical air mass at altitude z 
           msl : relative optical air mass at sea level
           height_corr : p/p0 */
        /* Outputs :
           deltaR : integral Rayleigh optical thickness per unit mass */
        /* The procedure "deltaR_Page" provides the integral Rayleigh optical 
           thickness i.e. the optical thickness of a pure Rayleigh scattering 
           atmosphere, per unit mass, along a specified path length. Returns 0 if OK, 
           1 otherwise. */
        {
            int ier;
            double a, b, c, d, e, m2, m3, m4, press_corr;

            ier = 0;
            a = 6.625928;
            b = 1.92969;
            c = -0.170073;
            d = 0.011517;
            e = -0.000285;
            if (msl <= 0.0)
                deltaR = 0.0;
            else
            {
                if (msl <= 20)
                {
                    m2 = msl * msl;
                    m3 = m2 * msl;
                    m4 = m3 * msl;

                    presscorr(height_corr, msl, out press_corr);
                    deltaR = (1.0 / (a + b * msl + c * m2 + d * m3 + e * m4)) / press_corr;
                }
                else  /* modification by Page (1996) m=20 gamma=3 degrees */
                    deltaR = 1.0 / (10.4 + 0.718 * m);
            }

            return (ier);
        }

        public static int presscorr(double height_corr, double msl, out double corr_press)

        /* Source :
           Page, J., 2001. (modified Rigollier et al. 2000 for pressure correction) */
        /* Inputs :
           height_corr : p/p0
           msl : relative optical air mass at sea level */
        /* Outputs :
           corr_press : pressure correction for deltaR, integral Rayleigh optical thickness per unit mass */
        /* Returns 0 if OK, 1 otherwise. */
        {
            int ier;
            double corr75, corr50, x, y;

            ier = 0;
            if (height_corr >= 0.99) { corr_press = 1.0; return ier; }

            corr75 = 1.248174 - 0.011997 * msl + 0.00037 * msl * msl;
            corr50 = 1.68219 - 0.03059 * msl + 0.00089 * msl * msl;

            /*  if(height_corr < 0.51) {*corr_press=corr50 ; return ier;} */

            if (height_corr > 0.75)
            {
                x = Math.Pow((height_corr - 0.75), 2.0);
                y = Math.Pow((height_corr - 1.0), 2.0);
                corr_press = (x * 1.0 + y * corr75) / (x + y);
                return ier;
            }
            x = Math.Pow((height_corr - 0.75), 2.0);
            y = Math.Pow((height_corr - 0.50), 2.0);
            corr_press = (x * corr50 + y * corr75) / (x + y);
            return (ier);
        }

        /******************************/
        /* LINKE TURBIDITY FACTOR, TL */
        /******************************/

        public static int TL_AM2_Angstrom(int month_number, double Angstrom_sum, out double TLAM2)

        /* Sources :
           Angstrom, 1924.
           Page, J., 1986.
           Solar Energy R&D in the European Community, Series F, Solar radiation data,
           Vol.3 D. Reidel Publishing Compagny, Dordrecht. 
           ESRA, 1984. Vol I. */
        /* Inputs :
           month_number : month number (1..12)
           Angstrom_sum : Angstrom sum (a+b) */
        /* Outputs :
           TLAM2 : Linke turbidity factor for air mass 2 */
        /* The procedure "TL_AM2_Angstrom" provides the Linke turbidity factor, 
           estimated at air mass 2, by way of a linear regression with the Angstrom 
           sum (a+b), i.e. the sum of the Angstrom regression coefficients. Returns 0 
           if OK, 1 otherwise. */
        {
            /* fm : factor chosen, according to the month of interest and to ... */
            /* ... monthly mean radiation calculations i.e. monthly averages values of 
               factor fm (CEC 1984 parametrization for monthly average clear skies) */
            double[] fm_average = new double[12] { 10.7, 13.2, 15.4, 17.1, 18.5, 16.9, 18.2, 17.0, 14.9, 12.9, 11.3, 9.5 };
            /* ... monthly mean maximum radiation predictions i.e. monthly maximum values 
               of factor fm (CEC 1984 parametrization for clearest skies within a month)*/
            /* static double fm_maximum[12]={13.5,13.9,13.4,15.4,15.4,15.4,14.9,16.3,14.6,13.2,12.8,11.6}; */
            int ier;

            ier = 1;
            TLAM2 = Double.NaN;
            if ((month_number > 0) && (month_number < 13))
            {
                ier = 0;
                TLAM2 = fm_average[month_number - 1] * (1.0 - Angstrom_sum);
                /* *TLAM2 = fm_maximum[month_number-1] * (1.0 - Angstrom_sum); */
            }

            return (ier);
        }

        public static int TL_AM2_ATI(double phi, double ATI, out double TLAM2)

        /* Sources :
           Page, J., 1986.
           Dogniaux and Lemoine, 1983.*/
        /* Inputs :
           phi : latitude (in radians, positive to North)
           ATI : Atmospheric Turbidity Index (estimated with the help of a table) */
        /* Outputs :
           TLAM2 : Linke turbidity factor for air mass 2 */
        /* The procedure "TL_AM2_ATI" provides the Linke turbidity factor, estimated at
           air mass 2, by using the Dogniaux and Lemoine parametrization from the 
           Atmospheric Turbidity Index. A subjective assessment of the type of 
           atmospheric conditions prevailing in the region of interest must be made 
           (table of ATI). Returns 0 if OK, 1 otherwise. */
        {
            const double deg_rad = (Math.PI / 180.0); /* converts decimal degrees into radians*/
            int ier;

            ier = 1;
            TLAM2 = Double.NaN;

            if (Math.Abs(phi) <= (Math.PI / 2.0))
            {
                ier = 0;
                TLAM2 = 22.76 + 3.071 * phi - 27.78 * ATI;
                /* Pb with the constant: 27.78 or 27.28 cf. ESRA2_06 */
            }

            return (ier);
        }

        public static int TLinke(double TLAM2, double gamma, out double TL)

        /* Sources :
           WMO, 1981.
           Meteorological aspects of the utilization of solar radiation as an energy 
           source. WMO-No. 557, Technical Note No. 172, Geneva, Switzerland, 298 pp.
           Page, J., 1986. */
        /* Inputs :
           TLAM2  : Linke turbidity factor for air mass 2
           gamma  : solar elevation angle (in radians) */
        /* Outputs :
           TL : Linke turbidity factor at a given solar altitude */
        /* The procedure "TLinke" estimates the Linke turbidity factor for any air 
           mass, TL(m), from its value at air mass 2, TL(2)=TLAM2, corrected for solar 
           elevation. Returns 0 if OK, 1 otherwise.
           NB : THIS PROCEDURE HAS TO BE USED IF THE PROCEDURE "deltaR_Kasten_old" IS 
           USED INSTEAD OF "deltaR_Page" */
        {
            const double deg_rad = (Math.PI / 180.0); /* converts decimal degrees into radians*/
            int ier;
            double sing = 0, a;

            ier = 1;
            TL = Double.NaN;
            if ((gamma >= 0.0) && (gamma <= (Math.PI / 2.0)))
            {
                ier = 0;
                sing = Math.Sin(gamma);
            }

            if ((gamma <= 0.0) || (TLAM2 < 1.0))
                TL = 0.0;
            else
            {
                a = (0.85 - 2.25 * sing + 1.11 * sing * sing);
                if (TLAM2 >= 2.5)
                    TL = TLAM2 - a;
                else
                    TL = TLAM2 - a * (TLAM2 - 1.0) / 1.5;

                /* To verify */
                /*
                a = TLAM2 - (0.85 - 2.25 * sing + 1.11 * sing * sing);
                if(TLAM2 >= 2.5)
                    *TL = a;
                else
                    *TL = a * (TLAM2 - 1.0) / 1.5;
                */

            }

            return (ier);
        }

        /*********************************************/
        /*                                           */
        /* C L E A R  S K Y  M O D E L  O F  P A G E */
        /*                                           */
        /*********************************************/

        public static double Beam_trans(double m, double deltaR, double TLAM2)

        /* Inputs :
           m      : relative optical air mass
           deltaR : integral Rayleigh optical thickness
           TLAM2  : Linke turbidity factor for air mass 2 */
        /* return : clear sky beam transmittance for beam radiation
                    if m = height_corr, beam transmittance at zenith Trb */
        {
            return Math.Exp(-0.8662 * TLAM2 * m * deltaR);
        }

        public static double Diff_trans(double TLAM2_corr)

        /* Inputs :
           TLAM2_corr  : Linke turbidity factor for air mass 2 corrected for p/p0 */
        /* return : clear sky transmittance at zenith for diffuse radiation, Trd */
        {
            return -1.5843e-02 + 3.0543e-02 * TLAM2_corr + 3.797e-04 * Math.Pow(TLAM2_corr, 2.0);
        }

        public static double Diff_ang_func(double TLAM2_corr, double gamma, double Trd)

        /* Inputs :
           TLAM2_corr  : Linke turbidity factor for air mass 2 corrected for p/p0 
           gamma  : solar elevation angle (in radians)
           Trd    : transmittance at zenith for diffuse radiation */
        /* return : diffuse angular function, Fd */
        {
            double sing, TLAM2_2, A0, A1, A2;
            double[][] a = new double[3][]{ new double[3] { 2.6463e-01,-6.1581e-02, 3.1408e-03},
                                  new double[3] { 2.0402    , 1.8945e-02,-1.1161e-02},
                                  new double[3] {-1.3025    , 3.9231e-02, 8.5079e-03} };

            TLAM2_2 = TLAM2_corr * TLAM2_corr;
            A0 = a[0][0] + a[0][1] * TLAM2_corr + a[0][2] * TLAM2_2;
            A1 = a[1][0] + a[1][1] * TLAM2_corr + a[1][2] * TLAM2_2;
            A2 = a[2][0] + a[2][1] * TLAM2_corr + a[2][2] * TLAM2_2;
            if ((Trd * A0) < 2.0e-03)
                A0 = 2.0e-03 / Trd;
            sing = Math.Sin(gamma);
            return A0 + A1 * sing + A2 * sing * sing;
        }

        public static int Diff_ang_func_irradiation(double phi, double delta, double omega1, double omega2, double TLAM2_corr, double Trd, out double Fd)

        /* Inputs :
           phi    : latitude of site (in radians, positive to North)
           delta  : solar declination angle (in radians)
           omega1 : solar hour angle at beginning of the time period (in radians)
           omega2 : solar hour angle at end of the time period (in radians) 
           TLAM2_corr : Linke turbidity factor for air mass 2
           Trd    : transmittance at zenith for diffuse radiation */
        /* Outputs :
           Fd      : diffuse angular function for irradiation */
        /* Return  : 0 if no error, 1 otherwise */
        {
            int ier = 0;
            Fd = Double.NaN;
            double TLAM2_2, omega_sr, omega_ss, A0, A1, A2, D0, D1, D2;
            double cospd, sinpd, cos2pd, sin2pd;
            double[][] a = new double[3][]{ new double[3] { 2.6463e-01,-6.1581e-02, 3.1408e-03},
                        new double[3]{ 2.0402    , 1.8945e-02,-1.1161e-02},
                        new double[3]{-1.3025    , 3.9231e-02, 8.5079e-03} };
            TLAM2_2 = TLAM2_corr * TLAM2_corr;
            cospd = Math.Cos(phi) * Math.Cos(delta);
            sinpd = Math.Sin(phi) * Math.Sin(delta);
            cos2pd = cospd * cospd;
            sin2pd = sinpd * sinpd;
            A0 = a[0][0] + a[0][1] * TLAM2_corr + a[0][2] * TLAM2_2;
            A1 = a[1][0] + a[1][1] * TLAM2_corr + a[1][2] * TLAM2_2;
            A2 = a[2][0] + a[2][1] * TLAM2_corr + a[2][2] * TLAM2_2;
            if ((Trd * A0) < 2.0e-03)
                A0 = 2.0e-03 / Trd;
            D0 = A0 + A1 * sinpd + A2 * sin2pd + 0.5 * A2 * cos2pd;
            D1 = A1 * cospd + 2.0 * A2 * sinpd * cospd;
            D2 = 0.25 * A2 * cos2pd;
            ier = SolarCalculator.sunrise_hour_angle(phi, delta, 0.0, out omega_sr, out omega_ss);
            if (ier != 0) return (ier);

            if ((omega2 < (omega_sr - (Math.PI / 12.0))) || (omega1 > (omega_ss + (Math.PI / 12.0))))
                Fd = 0.0;
            else
                Fd = D0 * (omega2 - omega1) + D1 * (Math.Sin(omega2) - Math.Sin(omega1)) +
                 D2 * (Math.Sin(2.0 * omega2) - Math.Sin(2.0 * omega1));
            return ier;
        }

        public static double Beam_ang_func(double phi, double delta, double omega_SR, double omega_SS, double omega1, double omega2, double TLAM2, double height_corr)

        /* Inputs :
           phi      : latitude of site (in radians, positive to North)
           delta    : solar declination angle (in radians)
           omega_SR : solar hour angle at sunrise (in radians)
           omega_SS : solar hour angle at sunset (in radians)
           omega1   : solar hour angle at beginning of the time period (in radians)
           omega2   : solar hour angle at end of the time period (in radians) 
           TLAM2    : Linke turbidity factor for air mass 2
           height_corr : p/p0 
        /* return Fb : Beam angular function for beam irradiation */
        {
            /* Coefficients for the beam angular function, Fb */
            double[][][] L = new double[3][][]{new double[3][]{new double[3] {-1.7349e-02,-5.8985e-03, 6.8868e-04},
                                                               new double[3] {-8.2193e-03, 4.5643e-04, 6.7916e-05},
                                                               new double[3]{-1.1656e-03, 1.8408e-04, -4.8754e-07} },
                                               new double[3][]{new double[3]{ 1.0258    ,-1.2196e-01, 1.9229e-03},
                                                               new double[3] { 8.9233e-01,-1.9991e-01, 9.9741e-03},
                                                               new double[3] { 7.4095e-01,-2.2427e-01, 1.5314e-02} },
                                               new double[3][]{new double[4]{-7.2178e-03, 1.3086e-01,-2.8405e-03, 0.0},
                                                               new double[4]  { 2.5428e-01, 2.6140e-01,-1.7020e-02, 0.0},
                                                               new double[4]  { 3.4959e-01, 7.2313e-01,-1.2305e-01, 5.9194e-03} } };

            double gamma_noon, cospd, sinpd, cos2pd, sin2pd;
            double Fb, b0, B1, B2, C0, C1, C2;
            int i = 0;

            gamma_noon = (Math.PI / 2.0) - phi + delta; /* solar elevation angle at noon */
            if (gamma_noon <= (Math.PI / 6.0)) i = 1;
            if (gamma_noon < (Math.PI / 12.0)) i = 2;

            C0 = L[0][i][0] + L[0][i][1] * TLAM2 * height_corr + L[0][i][2] * Math.Pow(TLAM2 * height_corr, 2.0);
            C1 = L[1][i][0] + L[1][i][1] * TLAM2 * height_corr + L[1][i][2] * Math.Pow(TLAM2 * height_corr, 2.0);
            C2 = L[2][i][0] + L[2][i][1] * TLAM2 * height_corr + L[2][i][2] * Math.Pow(TLAM2 * height_corr, 2.0) +
                L[2][i][3] * Math.Pow(TLAM2 * height_corr, 3.0);

            cospd = Math.Cos(phi) * Math.Cos(delta);
            sinpd = Math.Sin(phi) * Math.Sin(delta);
            cos2pd = cospd * cospd;
            sin2pd = sinpd * sinpd;

            b0 = C0 + C1 * sinpd + C2 * sin2pd + 0.5 * C2 * cos2pd;
            B1 = C1 * cospd + 2.0 * C2 * sinpd * cospd;
            B2 = 0.25 * C2 * cos2pd;

            if ((omega2 < omega_SR) || (omega1 > omega_SS))
                Fb = 0.0;
            else
            {
                if (omega1 < omega_SR)

                    omega1 = omega_SR;
                if (omega2 > omega_SS)
                    omega2 = omega_SS;
                Fb = b0 * (omega2 - omega1) + B1 * (Math.Sin(omega2) - Math.Sin(omega1)) +
              B2 * (Math.Sin(2.0 * omega2) - Math.Sin(2.0 * omega1));
            }
            return Fb;
        }

        public static int Bc_model5_direct(double eccentricity, double m, double deltaR, double TLAM2, out double Bc)

        /* Inputs :
           eccentricity : correction for Earth orbit eccentricity
           m            : relative optical air mass
           deltaR       : integral Rayleigh optical thickness
           TLAM2        : Linke turbidity factor for air mass 2 */
        /* Outputs :
           Bc : clear sky direct beam irradiance (in W/m2) */
        /* The procedure provides the cloudless sky direct beam irradiance (in W/m2)
           This procedure returns 0 if OK, 1 otherwise. */
        {
            int ier = 0;
            if ((TLAM2 <= 1.0) || (m <= 0.0))
                Bc = 0.0;
            else
                Bc = SolarCalculator.I0 * eccentricity * Beam_trans(m, deltaR, TLAM2);

            return (ier);
        }

        public static int Bc_model5_irradiance(double gamma, double eccentricity, double m, double deltaR, double TLAM2, out double Bc)

        /* Inputs :
           gamma        : solar elevation angle (in radians)
           eccentricity : correction for Earth orbit eccentricity
           m            : relative optical air mass
           deltaR       : integral Rayleigh optical thickness
           TLAM2        : Linke turbidity factor for air mass 2 */
        /* Outputs :
           Bc : clear sky beam horizontal irradiance (in W/m2) */
        /* The procedure provides the cloudless sky beam horizontal irradiance (in W/m2)
           This procedure returns 0 if OK, 1 otherwise. */
        {
            int ier;
            ier = Bc_model5_direct(eccentricity, m, deltaR, TLAM2, out Bc);
            if (ier == 0 && Bc > 0.0)
            {
                if (gamma > 0.0) Bc = Bc * Math.Sin(gamma); else Bc = 0.0;
            }
            return (ier);
        }

        public static int Bc_model5(double z, double eccentricity, double TLAM2, double phi, double delta, double omega_SR, double omega_SS, double omega1, double omega2, out double Bc)

        /* Inputs :
           z            : station height above sea level (in meters)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           delta        : solar declination angle (in radians)
           omega1       : solar hour angle at beginning of the time period (in radians)
           omega2       : solar hour angle at end of the time period (in radians) */
        /* Outputs :
           Bc : clear sky beam horizontal irradiation (in Wh/m2) */
        /* The procedure provides the cloudless sky beam horizontal irradiation (in Wh/m2)
           This procedure returns 0 if OK, 1 otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier;
            double dt, height_corr, deltaR, Trb;

            ier = 1;
            Bc = 0.0;
            if ((TLAM2 < 1.0) || (omega2 <= omega1))
            {
                ier = 0;
            }
            else
            {
                dt = SolarCalculator.Dl / (2.0 * Math.PI);
                station_height_correction(z, out height_corr);
                ier = deltaR_Page((double)1.0, height_corr, height_corr, out deltaR);
                if (ier != 0) return (ier);
                /* gamma = 90 deg i.e. m = 1 at sea level or m = p/p0 at higher station 
               heights --> air_mass_Kasten((Math.PI/2.0),z,&m) */

                /* Transmittance function for beam radiation at zenith Trb */
                Trb = Beam_trans(height_corr, deltaR, TLAM2);

                /* Clear sky beam horizontal irradiation : Bc */
                Bc = dt * SolarCalculator.I0 * eccentricity * Trb * Beam_ang_func(phi, delta, omega_SR, omega_SS, omega1, omega2, TLAM2, height_corr);
                if (Bc < 0.0) Bc = 0.0;
            }
            return (ier);
        }

        public static int Dc_model5_irradiance(double z, double gamma, double eccentricity, double TLAM2, out double Dc)

        /* Inputs :
           gamma        : solar elevation angle (in radians)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2 */
        /* Outputs :
           Dc : clear sky diffuse horizontal irradiance (in W/m2) */
        /* The procedure provides the clear sky diffuse horizontal irradiance (in W/m2)
           This procedure returns 0 if OK, 1 otherwise. */
        {
            int ier;
            double sing, Trd, Fd, height_corr, TLAM2_corr;

            ier = 1;
            if ((TLAM2 < 1.0) || (gamma < -(Math.PI / 12.0)))
            {
                ier = 0;
                Dc = 0.0;
            }
            else
            {
                ier = 0;
                sing = Math.Sin(gamma);
                /* Transmittance at zenith for diffuse radiation : Trd */
                station_height_correction(z, out height_corr);
                TLAM2_corr = TLAM2 * height_corr;
                Trd = Diff_trans(TLAM2_corr);
                Fd = Diff_ang_func(TLAM2_corr, gamma, Trd);
                Dc = SolarCalculator.I0 * eccentricity * Trd * Fd;
                if (Dc < 0.0) Dc = 0.0;
            }
            return (ier);
        }

        public static int Dc_model5(double z, double eccentricity, double TLAM2, double phi, double delta, double omega1, double omega2, out double Dc)
        /* Inputs :
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           delta        : solar declination angle (in radians)
           omega1       : solar hour angle at beginning of the time period (in radians)
           omega2       : solar hour angle at end of the time period (in radians) */
        /* Outputs :
           Dc : clear sky diffuse horizontal irradiation (in Wh/m2) */
        /* The procedure provides the cloudless sky diffuse horizontal irradiation (in Wh/m2)
           This procedure returns 0 if OK, 1 otherwise. */
        {
            int ier = 1;
            double dt, Trd, Fd, height_corr, TLAM2_corr;

            if ((TLAM2 < 1.0) || (omega2 <= omega1))
            {
                ier = 0;
                Dc = 0.0;
            }
            else
            {
                ier = 0;
                dt = SolarCalculator.Dl / (2.0 * Math.PI);
                /* Transmittance at zenith for diffuse radiation : Trd */
                station_height_correction(z, out height_corr);
                TLAM2_corr = TLAM2 * height_corr;
                Trd = Diff_trans(TLAM2_corr);
                /* Diffuse angular function : Fd */
                ier = Diff_ang_func_irradiation(phi, delta, omega1, omega2, TLAM2_corr, Trd, out Fd);

                Dc = dt * SolarCalculator.I0 * eccentricity * Trd * Fd;
                if (Dc < 0.0) Dc = 0.0;
            }
            return (ier);
        }

        public static int Gc_model5_irradiance(double gamma, double eccentricity, double TLAM2, double z, out double Bc, out double Dc, out double Gc)

        /* Inputs :
           gamma        : solar elevation angle (in radians)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           z            : station height above sea level (in meters) */
        /* Outputs :
           Bc : clear sky beam horizontal irradiance (in W/m2)
           Dc : clear sky diffuse horizontal irradiance (in W/m2)
           Gc : clear sky global horizontal irradiance (in W/m2) */
        /* Model of Page. The procedure "Gc_model5_irradiance" provides the clear sky 
           global horizontal irradiance (in W/m2) which is the sum of the clear sky 
           beam horizontal irradiance (in W/m2) and of the clear sky diffuse horizontal
           irradiance (in W/m2). This procedure returns 0 if OK, 1 otherwise. */
        {
            int ier;
            double m, msl, deltaR, height_corr;
            msl = deltaR = Double.NaN;

            ier = 1;
            Bc = 0.0;
            Dc = 0.0;
            Gc = 0.0;
            if ((TLAM2 < 1.0) || (gamma < -(Math.PI / 12.0)))
            {
                ier = 0;
            }
            else
            {
                station_height_correction(z, out height_corr);
                ier = air_mass_Kasten(gamma, z, out m);
                if (ier == 0) ier = air_mass_Kasten(gamma, (double)0.0, out msl);/* modif ML 02/02 */
                if (ier == 0) ier = deltaR_Page(msl, m, height_corr, out deltaR);
                if (ier != 0) return (ier);
                ier = Bc_model5_irradiance(gamma, eccentricity, m, deltaR, TLAM2, out Bc);
                if (ier == 0) ier = Dc_model5_irradiance(z, gamma, eccentricity, TLAM2, out Dc); /* modif L. Wald Sep 2014 */

                /* CLEAR SKY GLOBAL HORIZONTAL IRRADIANCE */
                if (ier == 0) Gc = Bc + Dc;  /* modif L. Wald Sep 2014 */
            }

            return (ier);
        }

        public static int Gc_model5(double z, double eccentricity, double TLAM2, double phi, double delta, double omega_SR, double omega_SS, double omega1, double omega2, out double Bc, out double Dc, out double Gc)

        /* Inputs :
           z            : station height above sea level (in meters)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           delta        : solar declination angle (in radians)
           omega_SR     : sunrise solar hour angle (in radians)
           omega_SS     : sunset solar hour angle (in radians)
           omega1       : solar hour angle at beginning of the time period (in radians)
           omega2       : solar hour angle at end of the time period (in radians) */
        /* Outputs :
           Bc : clear sky beam horizontal irradiation (in Wh/m2)
           Dc : clear sky diffuse horizontal irradiation (in Wh/m2)
           Gc : clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gc_model5" provides the clear sky global 
           horizontal irradiation (in Wh/m2) which is the sum of the clear sky beam 
           horizontal irradiation (in Wh/m2) and of the clear sky diffuse horizontal 
           irradiation (in Wh/m2). This procedure returns 0 if OK, 1 otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier = 1;
            Bc = 0.0;
            Dc = 0.0;
            Gc = 0.0;

            if ((TLAM2 < 1.0) || (omega2 <= omega1) || (omega_SS <= omega_SR))
            {
                ier = 0;
            }
            else
            {
                /* CLEAR SKY GLOBAL HORIZONTAL IRRADIATION : Gc */
                ier = Dc_model5(z, eccentricity, TLAM2, phi, delta, omega1, omega2, out Dc);
                if (ier == 0) ier = Bc_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, omega1, omega2, out Bc);
                if (ier == 0) Gc = Bc + Dc;
            }
            return (ier);
        }

        public static int Gcd_model5(double z, double eccentricity, double TLAM2, double phi, double delta, double omega_SR, double omega_SS, out double Bcd, out double Dcd, out double Gcd)

        /* Inputs :
           z            : station height above sea level (in meters)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           delta        : solar declination angle (in radians)
           omega_SR     : sunrise solar hour angle (in radians)
           omega_SS     : sunset solar hour angle (in radians) */
        /* Outputs :
           Bcd : daily clear sky beam horizontal irradiation (in Wh/m2)
           Dcd : daily clear sky diffuse horizontal irradiation (in Wh/m2)
           Gcd : daily clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gcd_model5" provides the daily clear sky 
           global horizontal irradiation (in Wh/m2) which is the sum of the daily clear
           sky beam horizontal irradiation (in Wh/m2) and of the daily clear sky 
           diffuse horizontal irradiation (in Wh/m2). This procedure returns 0 if OK, 1
           otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier = 1;
            Bcd = 0.0;
            Dcd = 0.0;
            Gcd = 0.0;
            if ((TLAM2 < 1.0) || (omega_SS <= omega_SR))
            {
                ier = 0;
            }
            else
            {
                /* CLEAR SKY GLOBAL HORIZONTAL IRRADIATION : Gc */
                ier = Dc_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, out Dcd);
                if (ier == 0) ier = Bc_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, omega_SR, omega_SS, out Bcd);
                if (ier == 0) Gcd = Bcd + Dcd;
            }
            return (ier);
        }

        public static int Gch24_model5(double z, double eccentricity, double TLAM2, double phi, double delta, double omega_SR, double omega_SS, out double[] Bch, out double[] Dch, out double[] Gch)

        /* Inputs :
           z            : station height above sea level (in meters)
           eccentricity : correction for Earth orbit eccentricity
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           delta        : solar declination angle (in radians)
           omega_SR     : sunrise solar hour angle (in radians)
           omega_SS     : sunset solar hour angle (in radians) */
        /* Outputs :
           Bch[1..24] : 24 hourly clear sky beam horizontal irradiation (in Wh/m2)
           Dch[1..24] : 24 hourly clear sky diffuse horizontal irradiation (in Wh/m2)
           Gch[1..24] : 24 hourly clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gch24_model5" provides the 24 hourly clear sky
           global horizontal irradiation (in Wh/m2) as well as the 24 hourly clear sky 
           beam horizontal irradiation (in Wh/m2) and the 24 hourly clear sky diffuse 
           horizontal irradiation (in Wh/m2). This procedure returns 0 if OK, 1 
           otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier, i;
            double t1, t2, omega1, omega2;
            Bch = new double[24];
            Dch = new double[24];
            Gch = new double[24];

            ier = 1;

            if ((TLAM2 < 1.0) || (omega_SS <= omega_SR))
            {
                ier = 0;
            }
            else
            {
                for (i = 0; i < 24; i++)
                {
                    t1 = ((i + 1) - 1.0);
                    ier = SolarCalculator.solar_hour_angle(t1, out omega1);
                    if (ier != 0) return (ier);
                    t2 = (i + 1);
                    ier = SolarCalculator.solar_hour_angle(t2, out omega2);
                    if (ier != 0) return (ier);

                    /* 24 HOURLY CLEAR SKY GLOBAL HORIZONTAL IRRADIATION, Gch[1..24] */
                    ier = Gc_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, omega1, omega2, out Bch[i], out Dch[i], out Gch[i]);
                }
            }
            return (ier);
        }


        /************************************/
        /* MODEL OF PAGE : MONTHLY AVERAGES */
        /************************************/
        public static int Gcdm_model5(int month_number, int year_number, double z, double TLAM2, double phi, double lambda, double gamma_riset,
                out double[] Bch_m, out double[] Dch_m, out double[] Gch_m, out double Bcd_m, out double Dcd_m, out double Gcd_m)


        /* Inputs :
           month_number : month number (1..12)
           year_number  : year number (4 digits)
           z            : station height above sea level (in meters)
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           lambda       : longitude of site (in radians, positive to East)
           gamma_riset  : solar elevation near sunrise/sunset:
                          - set to  0.0 for astronomical sunrise/sunset
                  - set to -1.0 for refraction corrected sunrise/sunset. */
        /* Outputs :      monthly average of...
           Bch_m[1..24] : 24 hourly clear sky beam horizontal irradiation (in Wh/m2)
           Dch_m[1..24] : 24 hourly clear sky diffuse horizontal irradiation (in Wh/m2)
           Gch_m[1..24] : 24 hourly clear sky global horizontal irradiation (in Wh/m2) 
           Bcd_m        : daily clear sky beam horizontal irradiation (in Wh/m2)
           Dcd_m        : daily clear sky diffuse horizontal irradiation (in Wh/m2)
           Gcd_m        : daily clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gcdm_model5" computes directly the monthly 
           average values of the 24 hourly clear sky beam, diffuse and global 
           horizontal irradiation (in Wh/m2) as well as monthly average values of the 
           daily clear sky beam, diffuse and global horizontal irradiation (in Wh/m2).
           This procedure returns 0 if OK, 1 otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier, i, day_of_month, number_days_month, julian_day;
            double day_angle, eccentricity, delta, omega_SR, omega_SS;
            double Bcd, Dcd, Gcd;

            double[] Bch = new double[24];
            double[] Dch = new double[24];
            double[] Gch = new double[24];

            Bch_m = new double[24];
            Dch_m = new double[24];
            Gch_m = new double[24];

            double nbd_m;

            /* Initialization */
            Bcd_m = 0.0;
            Dcd_m = 0.0;
            Gcd_m = 0.0;


            ier = 1;
            day_angle = delta = eccentricity = omega_SR = omega_SS = Double.NaN;
            ier = SolarCalculator.nbdays_month(year_number, month_number, out number_days_month);
            if (ier != 0) return (ier);

            for (day_of_month = 1; day_of_month <= number_days_month; day_of_month++)
            {
                ier = SolarCalculator.make_julian_day(day_of_month, month_number, year_number, out julian_day);
                if (ier == 0)
                    ier = SolarCalculator.Day_Angle(julian_day, out day_angle);
                if (ier == 0)
                    ier = SolarCalculator.corr_distance(day_angle, out eccentricity);
                if (ier == 0)
                    ier = SolarCalculator.declination_sun(year_number, julian_day, lambda, out delta);

                /* Computation of sunrise and sunset hour angles */
                if (ier == 0)
                    ier = SolarCalculator.sunrise_hour_angle(phi, delta, gamma_riset, out omega_SR, out omega_SS);

                /* OR OTHER SOLUTION */
                /* ier = solar_parameters_day(day_of_month,month_number,year_number,phi,lambda,gamma_riset,&day_angle,&delta,&omega_ss,&S0,&eccentricity,&G0d,G0h);
                if(ier != 0) return(ier); */

                if (ier == 0)
                    ier = Gch24_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, out Bch, out Dch, out Gch);
                if (ier != 0) return (ier);
                ier = Gcd_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, out Bcd, out Dcd, out Gcd);
                if (ier != 0) return (ier);

                for (i = 0; i < 24; i++)
                {
                    Bch_m[i] = Bch_m[i] + Bch[i];
                    Dch_m[i] = Dch_m[i] + Dch[i];
                    Gch_m[i] = Gch_m[i] + Gch[i];
                }
                Bcd_m = Bcd_m + Bcd;
                Dcd_m = Dcd_m + Dcd;
                Gcd_m = Gcd_m + Gcd;
            }

            nbd_m = (double)number_days_month;
            for (i = 0; i < 24; i++)
            {
                Bch_m[i] = Bch_m[i] / nbd_m;
                Dch_m[i] = Dch_m[i] / nbd_m;
                Gch_m[i] = Gch_m[i] / nbd_m;
            }
            Bcd_m = Bcd_m / nbd_m;
            Dcd_m = Dcd_m / nbd_m;
            Gcd_m = Gcd_m / nbd_m;

            return (ier);
        }

        /* BY USING A PONDERATION WITH ONLY THE POSITIVE VALUES OF IRRADIATION 
           COMPONENTS */
        public static int Gcdm_model5_2(int month_number, int year_number, double z, double TLAM2, double phi, double lambda, double gamma_riset,
                out double[] Bch_m, out double[] Dch_m, out double[] Gch_m, out double Bcd_m, out double Dcd_m, out double Gcd_m)

        /* Inputs :
           month_number : month number (1..12)
           year_number  : year number (4 digits)
           z            : station height above sea level (in meters)
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           lambda       : longitude of site (in radians, positive to East)
           gamma_riset  : solar elevation near sunrise/sunset:
                          - set to  0.0 for astronomical sunrise/sunset
                  - set to -1.0 for refraction corrected sunrise/sunset. */
        /* Outputs :      monthly average of...
           Bch_m[1..24] : 24 hourly clear sky beam horizontal irradiation (in Wh/m2)
           Dch_m[1..24] : 24 hourly clear sky diffuse horizontal irradiation (in Wh/m2)
           Gch_m[1..24] : 24 hourly clear sky global horizontal irradiation (in Wh/m2) 
           Bcd_m        : daily clear sky beam horizontal irradiation (in Wh/m2)
           Dcd_m        : daily clear sky diffuse horizontal irradiation (in Wh/m2)
           Gcd_m        : daily clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gcdm_model5" computes directly the monthly 
           average values of the 24 hourly clear sky beam, diffuse and global 
           horizontal irradiation (in Wh/m2) as well as monthly average values of the 
           daily clear sky beam, diffuse and global horizontal irradiation (in Wh/m2).
           This procedure returns 0 if OK, 1 otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier, i, day_of_month, number_days_month, julian_day;
            double day_angle, eccentricity, delta, omega_SR, omega_SS;
            double Bcd, Dcd, Gcd;
            double nbd_m;
            double n_Bcd_m, n_Dcd_m, n_Gcd_m;

            double[] Bch = new double[24];
            double[] Dch = new double[24];
            double[] Gch = new double[24];

            double[] n_Bch_m = new double[24];
            double[] n_Dch_m = new double[24];
            double[] n_Gch_m = new double[24];

            Bch_m = new double[24];
            Dch_m = new double[24];
            Gch_m = new double[24];

            /* Initialization */
            Bcd_m = 0.0;
            Dcd_m = 0.0;
            Gcd_m = 0.0;


            n_Bcd_m = 0.0;
            n_Dcd_m = 0.0;
            n_Gcd_m = 0.0;


            ier = 1;
            day_angle = delta = eccentricity = omega_SR = omega_SS = Double.NaN;

            ier = SolarCalculator.nbdays_month(year_number, month_number, out number_days_month);
            if (ier != 0) return (ier);

            for (day_of_month = 1; day_of_month <= number_days_month; day_of_month++)
            {
                ier = SolarCalculator.make_julian_day(day_of_month, month_number, year_number, out julian_day);
                if (ier == 0)
                    ier = SolarCalculator.Day_Angle(julian_day, out day_angle);
                if (ier == 0)
                    ier = SolarCalculator.corr_distance(day_angle, out eccentricity);
                if (ier == 0)
                    ier = SolarCalculator.declination_sun(year_number, julian_day, lambda, out delta);

                /* Computation of sunrise and sunset hour angles */
                if (ier == 0)
                    ier = SolarCalculator.sunrise_hour_angle(phi, delta, gamma_riset, out omega_SR, out omega_SS);

                /* OR OTHER SOLUTION */
                /* ier = solar_parameters_day(day_of_month,month_number,year_number,phi,lambda,gamma_riset,&day_angle,&delta,&omega_ss,&S0,&eccentricity,&G0d,G0h);
                if(ier != 0) return(ier); */

                if (ier == 0)
                    ier = Gch24_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, out Bch, out Dch, out Gch);
                if (ier != 0) return (ier);
                ier = Gcd_model5(z, eccentricity, TLAM2, phi, delta, omega_SR, omega_SS, out Bcd, out Dcd, out Gcd);
                if (ier != 0) return (ier);

                for (i = 0; i < 24; i++)
                {
                    Bch_m[i] = Bch_m[i] + Bch[i];
                    Dch_m[i] = Dch_m[i] + Dch[i];
                    Gch_m[i] = Gch_m[i] + Gch[i];
                }
                Bcd_m = Bcd_m + Bcd;
                Dcd_m = Dcd_m + Dcd;
                Gcd_m = Gcd_m + Gcd;

                for (i = 0; i < 24; i++)
                {
                    if (Bch[i] > 0.0) n_Bch_m[i] = n_Bch_m[i] + 1.0;
                    if (Dch[i] > 0.0) n_Dch_m[i] = n_Dch_m[i] + 1.0;
                    if (Gch[i] > 0.0) n_Gch_m[i] = n_Gch_m[i] + 1.0;
                }
                if (Bcd > 0.0) n_Bcd_m = n_Bcd_m + 1.0;
                if (Dcd > 0.0) n_Dcd_m = n_Dcd_m + 1.0;
                if (Gcd > 0.0) n_Gcd_m = n_Gcd_m + 1.0;
            }

            nbd_m = (double)number_days_month;
            for (i = 0; i < 24; i++)
            {
                if (n_Bch_m[i] == 0.0)
                    Bch_m[i] = 0.0;
                else
                    Bch_m[i] = Bch_m[i] / n_Bch_m[i];
                if (n_Dch_m[i] == 0.0)
                    Dch_m[i] = 0.0;
                else
                    Dch_m[i] = Dch_m[i] / n_Dch_m[i];
                if (n_Gch_m[i] == 0.0)
                    Gch_m[i] = 0.0;
                else
                    Gch_m[i] = Gch_m[i] / n_Gch_m[i];
            }
            if (n_Bcd_m == 0.0)
                Bcd_m = 0.0;
            else
                Bcd_m = Bcd_m / n_Bcd_m;
            if (n_Dcd_m == 0.0)
                Dcd_m = 0.0;
            else
                Dcd_m = Dcd_m / n_Dcd_m;
            if (n_Gcd_m == 0.0)
                Gcd_m = 0.0;
            else
                Gcd_m = Gcd_m / n_Gcd_m;

            return (ier);
        }

        /********************************************************************/
        /* YEARLY AVERAGES OF MONTHLY SOLAR RADIATION COMPONENTS            */
        /* (LONG TERM MEANS OF MONTHLY MEANS OF HOURLY AND DAILY SUMS)      */
        /********************************************************************/
        public static int Gcdym_model5(int month_number, int year_start, int year_end, double z, double TLAM2, double phi, double lambda, double gamma_riset, out double[] Bch_ym, out double[] Dch_ym, out double[] Gch_ym, out double Bcd_ym, out double Dcd_ym, out double Gcd_ym)

        /* Inputs :
           month_number : month number (1..12)
           year_start   : starting year of the considered period (4 digits)
           year_end     : ending year of the considered period (4 digits)
           z            : station height above sea level (in meters)
           TLAM2        : Linke turbidity factor for air mass 2
           phi          : latitude of site (in radians, positive to North)
           lambda       : longitude of site (in radians, positive to East)
           gamma_riset  : solar elevation near sunrise/sunset:
                          - set to  0.0 for astronomical sunrise/sunset
                  - set to -1.0 for refraction corrected sunrise/sunset. */
        /* Outputs :      monthly average of...
           Bch_ym[1..24]: 24 hourly clear sky beam horizontal irradiation (in Wh/m2)
           Dch_ym[1..24]: 24 hourly clear sky diffuse horizontal irradiation (in Wh/m2)
           Gch_ym[1..24]: 24 hourly clear sky global horizontal irradiation (in Wh/m2) 
           Bcd_ym       : daily clear sky beam horizontal irradiation (in Wh/m2)
           Dcd_ym       : daily clear sky diffuse horizontal irradiation (in Wh/m2)
           Gcd_ym       : daily clear sky global horizontal irradiation (in Wh/m2) */
        /* Model of Page. The procedure "Gcdm_model5" computes the yearly averages of 
           monthly average values of the 24 hourly clear sky beam, diffuse and global 
           horizontal irradiation (in Wh/m2) as well as monthly average values of the 
           daily clear sky beam, diffuse and global horizontal irradiation (in Wh/m2).
           This procedure returns 0 if OK, 1 otherwise. */
        /* NOTE: To account for (possibly asymmetric) obstructions near the horizon, 
           the sunrise and sunset hour angles omega_SR and omega_SS are required 
           independently; the obstructions are supposed to hide only the beam radiation
           component while diffuse radiation continues to reach the site. */
        {
            int ier, i, year_number;
            double Bcd_m, Dcd_m, Gcd_m;
            double number_of_years;

            double[] Bch_m = new double[24];
            double[] Dch_m = new double[24];
            double[] Gch_m = new double[24];

            Bch_ym = new double[24];
            Dch_ym = new double[24];
            Gch_ym = new double[24];

            ier = 1;
            number_of_years = (double)(year_end - year_start + 1.0);

            /* Initialization */
            Bcd_ym = 0.0;
            Dcd_ym = 0.0;
            Gcd_ym = 0.0;


            for (year_number = year_start; year_number <= year_end; year_number++)
            {
                ier = Gcdm_model5(month_number, year_number, z, TLAM2, phi, lambda, gamma_riset, out Bch_m, out Dch_m, out Gch_m, out Bcd_m, out Dcd_m, out Gcd_m);
                if (ier != 0) return (ier);

                for (i = 0; i < 24; i++)
                {
                    Bch_ym[i] = Bch_ym[i] + Bch_m[i];
                    Dch_ym[i] = Dch_ym[i] + Dch_m[i];
                    Gch_ym[i] = Gch_ym[i] + Gch_m[i];
                }
                Bcd_ym = Bcd_ym + Bcd_m;
                Dcd_ym = Dcd_ym + Dcd_m;
                Gcd_ym = Gcd_ym + Gcd_m;
                /*
                printf("year_number= %4d  Gcd_m = %8.2f (Wh/m2)  Gcd_ym = %8.2f (Wh/m2)\n",year_number,Gcd_m,*Gcd_ym);
                for(i=0;i<24;i++) 
                  printf("hh= %2d (hours)\tGch_m = %8.2f (Wh/m2)\tGch_ym = %8.2f (Wh/m2)\n",i+1,Gch_m[i],Gch_ym[i]);
                */
            }

            for (i = 0; i < 24; i++)
            {
                Bch_ym[i] = Bch_ym[i] / number_of_years;
                Dch_ym[i] = Dch_ym[i] / number_of_years;
                Gch_ym[i] = Gch_ym[i] / number_of_years;
            }
            Bcd_ym = Bcd_ym / number_of_years;
            Dcd_ym = Dcd_ym / number_of_years;
            Gcd_ym = Gcd_ym / number_of_years;

            return (ier);
        }
    }
}


