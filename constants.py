from math import pi

#Physical constants
c = 3.e8 # Speed of light in m/s
free_space_imp = 120*pi # Impedence of free space
k = 1.38e-23 # Boltzmann's constant in J/K
density_water = 1e3 # Density of liquid water in kg m^-3

#Mass conversions
g_per_kg = 1000.

#Distance conversions
cm_per_m = 100.
mm_per_cm = 10.
mm_per_m = mm_per_cm * cm_per_m
m_per_km = 1000.
cm_per_in = 2.54
nmi_per_m = 1./1852.354944

#Angular conversion
deg_per_rad = 180./pi # Duh.

#Time conversions
sec_per_min = 60.
min_per_hour = 60.
sec_per_hour = sec_per_min * min_per_hour
hour_per_day = 24.
sec_per_day = sec_per_hour * hour_per_day

#Temperature conversions
degf_per_degc = 9./5.
kelvin_offset = 273.15
