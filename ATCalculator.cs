using System;

using OIECalculators;

namespace AT
{
    //based on Steadman 94 : http://www.bom.gov.au/jshess/docs/1994/steadman.pdf

    public class TAData
    {
        public double TIndoor { get; set; }
        public double TShade { get; set; }
        public double TSun {get;set;}
    }



    public static class ATCalculator
    {
        //DateTime t must be in UTC
        public static double GetSolarAltitudeInDegrees(double LatitudeInDegrees, double LongitudeInDegress, DateTime t)
        {
            double lonrad = LongitudeInDegress / 360 * 2 * Math.PI;
            double latrad = LatitudeInDegrees / 360 * 2 * Math.PI;

            int julian_date;
            SolarCalculator.make_julian_day(t.Day, t.Month, t.Year, out julian_date);

            double day_angle;//annual day angle
            SolarCalculator.Day_Angle(julian_date, out day_angle);

            double LAT; // sun local time
            SolarCalculator.UT_to_LAT(t.TimeOfDay.TotalHours, day_angle, lonrad, out LAT);

            double delta;//sun declination
            SolarCalculator.declination_sun(t.Year, julian_date, lonrad, out delta);

            double omega;//solar hour angle
            SolarCalculator.solar_hour_angle(LAT, out omega);

            double gamma, theta;//altitude angle and zenithal angle (complementary)
            SolarCalculator.elevation_zenith_sun(latrad, delta, omega, out gamma, out theta);

            return gamma * 360 / 2 / Math.PI;
        }

        public static double GetExcentricity(DateTime t)
        {
            int julian_date;
            SolarCalculator.make_julian_day(t.Day, t.Month, t.Year, out julian_date);

            double day_angle;//annual day angle
            SolarCalculator.Day_Angle(julian_date, out day_angle);

            double eccentricity;
            SolarCalculator.corr_distance(day_angle, out eccentricity);

            return eccentricity;
        }


        //returns solar angle in degrees just in case it's needed elsewhere
        public static double  GetESRAIrradiances(DateTime UTCTime, double LatitudeInDegrees, double LongitudeInDegrees, int AltitudeInMeters, double LinkeTurbidity, out double Beam, out double Diffuse, out double Total)
        {
   
            // ESRA model needs the Linke Turbidity factor of the atmosphere of the location, it can be found here :
            //http://www.soda-pro.com/fr/web-services/atmosphere/linke-turbidity-factor-ozone-water-vapor-and-angstroembeta

            //calculate solar altitude angle;
            double SunAltiInDegrees = GetSolarAltitudeInDegrees(LatitudeInDegrees, LongitudeInDegrees, UTCTime);
            double SunAltiInRadians = SunAltiInDegrees / 360 * 2 * Math.PI;

            //calculate earth orbit correction, earth is 7% closer to the sun in northern winter than in northern summer 
            double Eccentricity = GetExcentricity(UTCTime);

            //call ESRA model to get beam, diffuse and global illumination, given solar angle 
            ESRACalculator.Gc_model5_irradiance(SunAltiInRadians, Eccentricity, LinkeTurbidity, AltitudeInMeters, out Beam, out Diffuse, out Total);

            return SunAltiInDegrees;


        }

        //degrees in Celsius
        public static TAData Compute (double Temp, double DewPoint, double WindInMPerSec,DateTime UTCTime, double LatitudeInDegrees, double LongitudeInDegrees,int AltitudeInMeters,double LinkeTurbidity)
        {
            TAData res = new TAData();

            //Compute Vapor pressure
            //if only rh is avail you need to compute dew point first,
            double Pa = 6.11 * Math.Pow(10, 7.5 * DewPoint / (DewPoint + 237.3)) / 10; // /10 because Steadman formulas need kPa.

            // Indoor temp
            res.TIndoor=0.89*Temp + 3.82 * Pa - 2.56;     
            
            // Shade temp, we use the influence of wind
            res.TShade=Temp+3.3*Pa-0.70 * WindInMPerSec - 4;

            //now the fun part, we add clear sky illumination
            // we will compute beam and diffuse illuminations using ESRA model.

            // Beam illumination is QD in Steadman, Diffuse is Qd
            double QD, Qd, Qtotal;
            
            double SunAltiInDegrees =GetESRAIrradiances(UTCTime, LatitudeInDegrees, LongitudeInDegrees, AltitudeInMeters, LinkeTurbidity, out QD, out Qd, out Qtotal);
            double SunAltiInRadians = SunAltiInDegrees / 360 * 2 * Math.PI;

            //Calculate human projected area factor, it varies between 0.3 m2 when the sun is horizontal, to 0.1 m2 when it's vertical
            // using Steadman formula 
            double psi3 =     0.3014 
                            + 88E-5 * SunAltiInDegrees 
                            - 746E-7 * SunAltiInDegrees * SunAltiInDegrees 
                            + 435E-9 * SunAltiInDegrees * SunAltiInDegrees * SunAltiInDegrees;

            //now calculate Q1, direct radiation contribution
            double Q1;
            if (SunAltiInRadians > 0)
                Q1 = 0.5 * psi3 * QD / Math.Sin(SunAltiInRadians);
            else
                Q1 = 0;

            //Q2 is diffuse radiation contribution
            double Q2 = Qd / 7;

            //Q3 is radiation of the earth, using an albedo of 0.2
            double Q3 = Qtotal / 28;

            //now calculate the outgoing radiation from the body, considering a cloudless sky
            double AltitudeInKm = AltitudeInMeters / 1000.0;
            double Q4= (103 + Temp) * (1 - Math.Exp(-0.11 * AltitudeInKm) + Math.Exp(-0.11 * AltitudeInKm - 1.1)) * Math.Exp(-0.10 * Pa);

            // Total Qg
            double Qg = Q1 + Q2 + Q3 - Q4;


            //now we can compute apparent temperature
            res.TSun = Temp + 3.48 * Pa - 0.7 * WindInMPerSec + 0.7 * Qg / (WindInMPerSec + 10) - 4.25;


            //didn't do massive tests but Qg seems to be between -30 and +120, yielding a max +8K effect like described in the paper.

            return res;
        }

    }
}
