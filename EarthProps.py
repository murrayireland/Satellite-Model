# Earth settings

import math

Mass = 5.972186e24            # Mass [kg]
SemiMajorAxis = 6378137       # Equatorial radius [m]
Flattening = 1/298.257223563  # Flattening
AngularRate = 7292115e-11     # Angular velocity [rad/s] 
            
SemiMinorAxis = SemiMajorAxis*(1 - Flattening)
Eccentrcity = math.sqrt((SemiMajorAxis**2 - SemiMinorAxis**2)/SemiMajorAxis**2)
        