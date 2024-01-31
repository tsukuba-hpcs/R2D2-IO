# Input parameters
This file provides summary of input parameters and some suggestion for "ordinary " calculations.

## Basic parameters

- Mass:$1.001 M_\odot$
- Radius: $1.000 R_\odot$
- Luminosity: $1.005 L_\odot$
- Age: $4.600\times 10^{9} ~\mathrm{yr}$

## Photosphere

### Surface values

- Temperature: $5850.177 ~\mathrm{K}$
- Gravity: $2.740\times 10^{4} ~\mathrm{cm/s^2}$
- Density: $2.027\times 10^{-7} ~\mathrm{g/cm^3}$ 
- Pressure: $7.822\times 10^{4} ~\mathrm{dyn/cm^2}$ 
- Pressure scale height: $140.83 ~\mathrm{km}$
- Convection velocity (Maximum value):  3.48 km/s 
- Speed of sound: $7.95 ~\mathrm{km/s}$
- Convection time scale (at maximum $v_\mathrm{c}$): $3.44 ~\mathrm{min}$
- Equipartition magnetic field to 
  - Gas pressure: $1402.1 ~\mathrm{G} = \sqrt{8\pi p_\mathrm{surface}}$ 
  - Convection flow: $556.2 ~\mathrm{G} =  v\sqrt{4\pi\rho_\mathrm{surface}}$ 

### Photospheric calculation suggestion

#### Calculation domain

- Vertical (upper):  $R_* +  653.0 ~\mathrm{km}$ where  $\rho_\mathrm{surface}/\rho = 250$  
- Vertical (lower): $R_*  -4501.2 ~\mathrm{km}$ where   $\rho_\mathrm{surface}/\rho = 1/500$  
- Horizontal: $5633.3 ~\mathrm{km}~=4H_p\times 10$

#### Grid spacing (typical)
- Vetical: $35.2 ~\mathrm{km} =  H_p/4$ 
- Horizontal: $70.4 ~\mathrm{km} =  H_p/2$ 

#### Time
- Output cadence: $20.6 ~\mathrm{s} =  4H_p/v_\mathrm{c}/10$ 
- Calculation Duration: $51.5 ~\mathrm{min} =  4H_p/v_\mathrm{c}\times15$ 

### With sunspot
- Magnetic field strength = $3505.3 ~\mathrm{G} =  2.5\sqrt{8\pi p_\mathrm{surf}}$ 
- Sunspot radius = $11.3 ~\mathrm{Mm} = 80H_p$ 

## Deep convection zone

### Deep Convection zone values

- Location of base of the convection zone:  $R_\mathrm{bcz}$= 0.711 $R_\mathrm{*}$ 

Following values are at the middle of the convection zone $(R_* + R_\mathrm{bcz})/2$ 

- Temperature: $9.208\times 10^{5} ~\mathrm{K}$
- Gravity: $3.728\times 10^{4} ~\mathrm{cm/s^2}$
- Density:$5.100\times 10^{-2} ~\mathrm{g/cm^3}$ 
- Pressure scale height:  33.55 Mm 
- Convection velocity:  68.03 m/s 
- Convection time scale:  22.84 day 

### Deep CZ calculation suggestion

#### calculation domain

- Vertical (upper): $0.9558  R_*$  where   $\rho_\mathrm{bcz}/\rho = 30$  
- Vertical (lower): $0.7109 R_*=R_\mathrm{bcz}$  

- Horizontal domain should cover whole sphere

#### Grid spacing

- Vetical: $0.048 R_*=H_p/4$ 
- Horizongal grid spacing is arbitrary

#### Time
- Output cadence: $2.3 ~\mathrm{day} = 4H_p/v_\mathrm{c}/10$ 
- Calculation Duration: $342.6 ~\mathrm{day} = 4H_p/v_\mathrm{c}\times15$ 

## Appendix

- Solar mass: $1.988\times 10^{33} ~\mathrm{g}$
- Solar radius: $6.957\times 10^{10} ~\mathrm{cm}$
- Solar luminosity: $3.828\times 10^{33} ~\mathrm{erg/s}$
- Solar age: $4.570\times 10^{9} ~\mathrm{Gyr}$

