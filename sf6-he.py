import CoolProp.CoolProp as CP

# seen http://www.coolprop.org/dev/fluid_properties/Mixtures.html#mixtures


CP.apply_simple_mixing_rule('Helium', 'SF6', 'linear')

def properties_he_sf6(T,P,x):
	''' Fluid Property of He-SF6 mixture with x-He
	'''

	fluidstring = 'Helium[%f]&SF6[%f]'%(x,1-x)
	T = T + 273.15 # in K
	props = {}	
	props['rho']  = CP.PropsSI('Dmass','T',T,'P',P,fluidstring)
	props['lambda'] = CP.PropsSI('L','T',T,'P',P,fluidstring)
	props['cp']   = CP.PropsSI('CPMASS','T',T,'P',P,fluidstring)
	props['eta']  = CP.PropsSI('VISCOSITY','T',T,'P',P,fluidstring)
	props['alpha'] = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T,'P',P,fluidstring)
	props['nu'] = props['eta']/props['rho']
	props['kappa'] = props['lambda']/(props['cp']*props['rho'])
	props['Pr'] = props['nu']/props['kappa']
	return props 


T = 20 # temperature in degC
P = 18e5 # pressure in Pa

g  = 9.81 # gravity
H  = 4    # height in m
dT = 20   # dT in K

props = properties_he_sf6(T,P,0.5)

Ra = g*H**3*dT*props['alpha']/(props['nu']*props['kappa'])

print("Prandtl: %lf"%props['Pr'])
print("Rayleigh: %g"%Ra)

