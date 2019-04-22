import numpy as np
from scipy.integrate import quad

def f_S_lens(d, a, s):
	# all possible sites with radius [s] intersecting radius of accuracy [a] around a registered action and an observed point at distance [d] from this action
	# d: distance between observed point and point of registration of action
	# a: radius of accuracy of the action
	# s: radius of a site
	
	def fl1(d, r1, r2):
		# length of a circular segment created by intersecting two circles with distance d and radii r1 and r2
		# if r2 > sqrt(r1^2 - d^2)
		
		return 2*r2*np.arcsin(np.sqrt(-r2**4+(2*d**2+2*r1**2)*r2**2-d**4+2*r1**2*d**2-r1**4)/(2*d*r2))
	
	def fl2(d, r1, r2):
		# length of a circular segment created by intersecting two circles with distance d and radii r1 and r2
		# if r2 < sqrt(r1^2 - d^2)
		
		return 2*np.pi*r2-2*r2*np.arcsin(np.sqrt(-r2**4+(2*d**2+2*r1**2)*r2**2-d**4+2*r1**2*d**2-r1**4)/(2*d*r2))
	
	def fArea(d, r):
		# area of a lens created by intersecting two circles with distance d and equal radii r
		
		return 2*r**2*np.arccos(d/(2*r))-(d/2)*np.sqrt(4*r**2 - d**2)
	
	if d > 2*s + a:
		return 0
	
	elif d == 0:
		if 2*s <= a:
			return np.pi**2*s**4
		
		else:
			return (np.pi*(-2*a*s**2*np.sqrt(4*s**2-a**2)-a**3*np.sqrt(4*s**2-a**2)+8*np.arcsin(a/(2*s))*s**4+8*a**2*np.arccos(a/(2*s))*s**2))/4
	
	elif (d < a - 2*s):
		return np.pi**2*s**4
	
	elif d >= a:
		if d > 2*s - a:
			return quad(lambda R: fl1(d, a, R) * fArea(R, s), d - a, 2*s)[0]
		else:
			return quad(lambda R: fl1(d, a, R) * fArea(R, s),	d - a,	d + a)[0]
	
	elif d < a:
		I1 = I2 = I3 = 0
		lim2 = np.sqrt(a**2-d**2)
		if d > 2*s - a:
			if 2*s > lim2:
				I1 = quad(lambda R: fl1(d, a, R) * fArea(R, s), lim2, 2*s)[0]
				I2 = quad(lambda R: fl2(d, a, R) * fArea(R, s), a - d, lim2)[0]
				I3 = (-8*np.pi*np.arcsin((d-a)/(2*s))*s**4+np.sqrt(4*s**2-d**2+2*a*d-a**2)*((2*np.pi*d-2*np.pi*a)*s**2+np.pi*d**3-3*np.pi*a*d**2+3*np.pi*a**2*d-np.pi*a**3)+((-8*np.pi*d**2+16*np.pi*a*d-8*np.pi*a**2)*np.arccos((d-a)/(2*s))+8*np.pi**2*d**2-16*np.pi**2*a*d+8*np.pi**2*a**2)*s**2)/4
			elif 2*s > a - d:
				I2 = quad(lambda R: fl2(d, a, R) * fArea(R, s), a - d, 2*s)[0]
				I3 = (-8*np.pi*np.arcsin((d-a)/(2*s))*s**4+np.sqrt(4*s**2-d**2+2*a*d-a**2)*((2*np.pi*d-2*np.pi*a)*s**2+np.pi*d**3-3*np.pi*a*d**2+3*np.pi*a**2*d-np.pi*a**3)+((-8*np.pi*d**2+16*np.pi*a*d-8*np.pi*a**2)*np.arccos((d-a)/(2*s))+8*np.pi**2*d**2-16*np.pi**2*a*d+8*np.pi**2*a**2)*s**2)/4
			else:
				I3 = np.pi**2*s**4
	
		else:
			I1 = quad(lambda R: fl1(d, a, R) * fArea(R, s), lim2, a + d)[0]
			I2 = quad(lambda R: fl2(d, a, R) * fArea(R, s), a - d, lim2)[0]
			I3 = (-8*np.pi*np.arcsin((d-a)/(2*s))*s**4+np.sqrt(4*s**2-d**2+2*a*d-a**2)*((2*np.pi*d-2*np.pi*a)*s**2+np.pi*d**3-3*np.pi*a*d**2+3*np.pi*a**2*d-np.pi*a**3)+((-8*np.pi*d**2+16*np.pi*a*d-8*np.pi*a**2)*np.arccos((d-a)/(2*s))+8*np.pi**2*d**2-16*np.pi**2*a*d+8*np.pi**2*a**2)*s**2)/4
		
		return I1 + I2 + I3

def f_S(a, s):
	# all possible sites with radius [s] intersecting radius of accuracy [a] around a registered action
	# a: radius of accuracy of the action
	# s: radius of a site	
	
	return np.pi**2*s**2*a**2

def f_s_approx(d, a, s):
	# approximation of f_S_lens / f_S
	
	if d > 2*s + a:
		return 0
	
	if 2*s <= a:
		n = s**2/a**2
	else:
		n = (8*np.arcsin(a/(2*s))*s**4+8*a**2*np.arccos(a/(2*s))*s**2-2*a*s**2*np.sqrt(4*s**2-a**2)-a**3*np.sqrt(4*s**2-a**2))/(4*np.pi*a**2*s**2)
	
	if d == 0:
		return n
	
	elif (d < a - 2*s):
		return s**2/a**2
	
	div = 1
	if max(a, s) > 1:
		div = max(a, s)
		a /= div
		s /= div
		if 2*s <= a:
			n = s**2/a**2
		else:
			n = (8*np.arcsin(a/(2*s))*s**4+8*a**2*np.arccos(a/(2*s))*s**2-2*a*s**2*np.sqrt(4*s**2-a**2)-a**3*np.sqrt(4*s**2-a**2))/(4*np.pi*a**2*s**2)
	
	if s < (a*0.3252 - 0.0965):
		l = s*0.195317 +(1.0197*a -0.0243)
	elif s < a*0.6451:
		l = 1.0587*a
	else:
		l = s*(0.0783*a**2 +0.0248*a +0.9968)+(0.1977*a**2 +0.1817*a +0.0258)
	
	if a > 0.5:
		if s < (0.7179*a +0.0801):
			k = a / (0.6451*s)
		else:
			k = 2.030303
	elif s < 0.3:
		k = a / (0.6451*s)
	elif s < (3.3*a +0.07):
		k = 2.030303
	else:
		k = 1.686868
	
	return n - n*(1-np.exp(-((d/div)/l)**k))

def f_t_UPD(t, t_d, sigma_d, sigma_s):
	# temporal component of the EDE function when using UPD-dated evidence
	# the ratio of possible sites, registered by an action which existed at time t and all such possible sites
	# t: observed date (temporal coordinate) in years
	# t_d: mean dating of the action
	# sigma_d: temporal uncertainty of the dating of the action in years
	# sigma_s: expected half-life of a site in years

	# dt = absolute distance in time between t and t_d 
	dt = abs(t - t_d)
	
	if dt >= sigma_d + 2*sigma_s:
		return 0
	
	elif dt < sigma_d - 2*sigma_s:
		return sigma_s/sigma_d
	
	elif dt > 2*sigma_s - sigma_d:
		return ((dt-sigma_d)*abs(dt-sigma_d)-4*sigma_s*dt+4*sigma_s**2+4*sigma_d*sigma_s)/(8*sigma_d*sigma_s)
	
	elif dt <= 2*sigma_s - sigma_d:
		return ((dt-sigma_d)*abs(dt-sigma_d)-dt**2-2*sigma_d*dt+8*sigma_d*sigma_s-sigma_d**2)/(8*sigma_d*sigma_s)

def f_t_NPD(t, sigma_s, p_dist, ts):
	# temporal component of the EDE function when using NPD-dated evidence
	# the ratio of possible sites, registered by an action which existed at time t and all such possible sites
	# t: observed date (temporal coordinate) in years
	# sigma_s: expected half-life of a site in years
	# p_dist: probability distribution resulting from calibration of the radiocarbon date
	# ts: calendar dates corresponding with the values of p_dist
	
	f_sigma_s_t = 2*sigma_s - np.abs(ts - t)
	f_sigma_s_t[f_sigma_s_t < 0] = 0

	return (p_dist*f_sigma_s_t).sum() / (p_dist.sum()*2*sigma_s)
