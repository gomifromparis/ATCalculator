using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ATCalculator
{
    class Program
    {
        static void Main(string[] args)
        {
            double Latitude = 48;
            double Longitude = 20;
            int Altitude = 0; //in m

            double Temperature = 20; //°c
            double DewPoint = 14;//°c
            double Wind = 2; // in m/s

            double TurbidityFactor = 3.2; //from http://www.soda-pro.com/fr/web-services/atmosphere/linke-turbidity-factor-ozone-water-vapor-and-angstroembeta)
            DateTime t = new DateTime(2019, 06, 01, 18, 15, 00).ToUniversalTime(); // Time must be UTC

            AT.TAData dat = AT.ATCalculator.Compute(Temperature, DewPoint, Wind, t, Latitude, Longitude, Altitude, TurbidityFactor) ;

            Console.WriteLine("Nominal:\t" + Temperature.ToString("0.00"));
            Console.WriteLine("Indoor:\t\t" + dat.TIndoor.ToString("0.00"));
            Console.WriteLine("Shaded:\t\t" + dat.TShade.ToString("0.00"));
            Console.WriteLine("Clear Sky:\t" + dat.TSun.ToString("0.00"));

        }
    }
}
