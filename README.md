# ATCalculator
This c# project will compute indoor, shaded and clear sky apparent temperatures following the Steadman 1994 paper. 
Results seem to give a max radiation effect of +8K which is in accordance with the paper.

Sun direct and diffuse illumations are calculated using a C# european ESRA model
Sun position calculations and the ESRA calculations are C# ports C sources available on the  the Ecole des Mines OIE Lab (see http://www.oie.mines-paristech.fr/Accueil/)

References:

Steadman :

1994: http://www.bom.gov.au/jshess/docs/1994/steadman.pdf

1984 : https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450(1984)023%3C1674:AUSOAT%3E2.0.CO;2

1979 : https://wonder.cdc.gov/wonder/help/Climate/Steadman1979.PDF

ESRA : 
https://hal.archives-ouvertes.fr/hal-00361373/document

http://www.oie.mines-paristech.fr/Valorisation/Outils/Solar-Geometry/

http://www.oie.mines-paristech.fr/Valorisation/Outils/Clear-Sky-Library/

http://www.youla.eu/PV/DOCS/PV_DIM.doc


Linke turbidity factors needed for the ESRA model from any site can be downloaded from here :
http://www.soda-pro.com/fr/web-services/atmosphere/linke-turbidity-factor-ozone-water-vapor-and-angstroembeta

