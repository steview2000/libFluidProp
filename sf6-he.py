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

	# calculate separation ratio  according to Cross and Hohenberg, Rev. Mod. Phys., (1993), page 999
	x2 = x*1.001
	x1 = x*0.999
	
	fluidstring1 = 'Helium[%f]&SF6[%f]'%(x1,1-x1)
	fluidstring2 = 'Helium[%f]&SF6[%f]'%(x2,1-x2)
	rho1  = CP.PropsSI('Dmass','T',T,'P',P,fluidstring1)
	rho2  = CP.PropsSI('Dmass','T',T,'P',P,fluidstring2)
	beta = (rho2-rho1)/(props['rho'] * (x2-x1))
	c  = beta/props['alpha']
	
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

