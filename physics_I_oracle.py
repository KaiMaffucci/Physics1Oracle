# imported math functions
from math import sin, cos, asin, acos, degrees, radians, sqrt, cbrt

"""
TODO:
 - anything else from exams/hw

TODO: test it all
maybe an exist() function (returns true if it isn't <None>) instead of all these '== None's ?
could put each chapter/unit into a class which can store data, and add a UI, but that's only a potential future endevour
"""

# MASTER HELP FUNCTION

def oracle_help(category = None):

    if category == "units":

        print("UNIT CONVERSIONS:")
        print("\tg_to_kg(g)")
        print("\tkm_to_m(km)")
        print("\tcm_to_m(cm)")
        print("\tm_to_cm(m)")
        print("\tm_to_km(m)")
        print("\tAU_to_m(AU)")
        print("\tdays_to_s(days)\n")

    elif category == "comps":

        print("COMPONENT SPLITTERS:")
        print("\tx component(x, theta)")
        print("\ty component(y, theta)\n")
    
    elif category == "vecs":

        print("VECTOR FUNCTIONS:")
        print("\tmag(A)")
        print("\tscalar_by_vector(n, A)")
        print("\tdot(A, B)")
        print("\tangle_between_vectors(A, B)")
        print("\tcross(A, B)")
        print("\tcenter(A, B)\n")

    elif category == "forces":

        print("FORCES / NEWTON\'S 2ND LAW BASICS")
        print("\tforce_info(F, m, a)")
        print("\ttension_help()")
        print("\tfriction_help()\n")

    elif category == "kine":
        
        print("KINEMATICS (linear & angular)")
        print("\tkinematics_help()\n")

    elif category == "energy":

        print("ENERGY")
        print("\tenergy_help()\n")
    
    elif category == "mom":

        print("MOMENTUM:")
        print("\tp_info(p, m, v)")
        print("\tp_help()\n")

    elif category == "centripedal" or category == "cent":
        
        print("CENTRIPEDAL FORCE")
        print("\tf_cp_info(F_cp, m, v, r)")
        print("\tmin_vertical_v_cp(r)")
        print("\tf_cp_help()")

    elif category == "angular" or category == "ang":

        print("ANGULAR SHIT:")
        print("\tangular_help()\n")
    
    elif category == "space":

        print("SPACE AND SHIT:")
        print("\tgrav_force_info(F, U, m1, m2, r)")
        print("\tplanet_info(g_planet, M_planet, R_planet)")
        print("\tescape_info(v_escape_planet, M_planet, R_planet)")
        print("\tplanet_period_info\n")
    
    elif category == "harmonics":

        print("SIMPLE HARMONIC MOTION:")
        print("\tperiod_info(T, m, k)")
        print("\thooke_info(m, k, x, a = g)")
        print("\tspring_energy_info(k, x, m, v, E_con)")
        print("\tperiod_info_alt(T, L, a = g)\n")

    else:

        print("Categories:")
        print("\tUNIT CONVERSTIONS (\'units\')")
        print("\tCOMPONENT SPLITTERS (\'comps\')")
        print("\tVECTOR FUNCTIONS (\'vecs'\)")
        print("\tFORCES (\'forces\')")
        print("\tKINEMATICS, linear and angular (\'kine\')")
        print("\tENERGY (\'energy\')")
        print("\tMOMENTUM (\'mom\')")
        print("\tANGULAR MOTION (\'angular\')")
        print("\tSPACE (\'space\')")
        print("\tSIMPLE HARMONIC MOTION (\'harmonic\')\n")



# CONSTANTS

from math import pi             # pi
g = 9.8                         # acceleration due to gravity (m/s**2)
G = 6.67 * 10**(-11)            # G gravity constant
r_earth = 6.37 * 10**6          # radius of earth (m)
M_earth = 5.972 * 10**24        # mass of earth (kg)



# UNIT CONVERSIONS


# g to kg
def g_to_kg(grams):
    return grams / 1000

# km to m
def km_to_m(km):
    return km * 1000

# cm to m
def cm_to_m(cm):
    return cm / 100

# m to cm
def m_to_cm(m):
    return m * 100

# m to km
def m_to_km(m):
    return m / 1000

# astronomical units (AU) to meters
def AU_to_m(AU): 
    return AU * 1.50 * 10**11

# days to seconds
def days_to_s(days):
    return 3600 * 24 * days



# X & Y COMPONENTS

"""
x component
theta is in degrees
"""
def xcomp(x, theta):
    return x * cos(degrees(theta))


"""
y component
theta is in degrees
"""
def ycomp(y, theta):
    return y * sin(degrees(theta))



# FUNCTIONS OF VECTORS

# checks if two given vectors' dimensions match or not (helper function)
def dim_mismatch(A, B):
    if (len(A) != len(B)): return True

# gets the magnitude of a vector (as a tuple or list)
def mag(A):
    
    abs_A = 0
    for dim in A:
        abs_A += dim**2
    
    abs_A = sqrt(abs_A)
    return abs_A


"""
Multiplies vector by scalar, returning new vector
Parameters: scalar, vector (tuple/list)
"""
def scalar_by_vector(n, A):

    new_A = list(A)
    for dim in range(len(new_A)): new_A[dim] *= n
    return tuple(new_A)


"""
Dot product
Parameters: two vectors of the same dimensions, as tuples/lists
"""
def dot(A, B):

    # if dimensions mismatch, function fails
    if dim_mismatch(A, B): exit

    # if they do, calculate & return dot product
    product = 0
    for dim in range(len(A)):
        product += A[dim] * B[dim]
    
    return product


"""
Finds angle between two vectors, in degrees
Parameters: two vectors of the same dimensions, as tuples/lists
"""
def angle_between_vecs(A, B):
    return degrees( acos(dot(A, B) / (mag(A)*mag(B))) )


"""
Finds cross product of two vectors
Parameters: two three-dimensional vectors, as tuples/lists
Tip: if you want to use this for vectors without a component of a dimension,
just input that value as 0
"""
def cross(A, B):
    return (+(A[1]*B[2] - A[2]*B[1]),
            -(A[0]*B[2] - A[2]*B[0]),
            +(A[0]*B[1] - A[1]*B[0]))


"""
Used mainly for things like finding center of mass
Returns dot product of two vectors / the sum of all items in the 2nd (eg. radii and mass, for center of mass)
Parameters: two vectors, as tuples or lists
"""
def center(A, B):
    return dot(A, B) / sum(B)



# FORCES / NEWTON'S 2ND LAW

"""
Gets info relating to the basic F = ma formula, based on what's <None>
Parameters: force (N), mass (kg), acceleration (a)
"""
def force_info(F, m, a):

    # calc whatever's unknown
    if F == None: F = m * a
    if m == None: m = F / a
    if a == None: a = F / m

    # print results after calc
    print("Force (N):", F)
    print("Mass (kg):", m)
    print("Acceleration (m/s**2):", a)


# Explains how tension works
def tension_help():
    
    print("TENSION FORCES")
    print("Tension exerted is a sort of \"counter-force\" to the force of gravity (or whatever other force may be at play)")
    print("Tension goes in both directions, but those extra tension forces usually cancel each other out.")
    print("net force = force due to gravity - tension force")
    print("Net force can be made negative if it's direction is being changed, say in a pulley system with two weights.")
    print("Two of these formulas can be constructed in a pulley system,")
    print("which can be used to solve for something like acceleration.\n")


# Sometimes you don't have the value of a force calculated already, so I left this as just a help function
def friction_help():

    print("FRICTION FORCES")
    print("Friction force = (coefficient of friction)(normal force)")
    print("Normal force = y component of force due to gravity (usually)")
    print("Trying to get an object moving? That force has to be greater than the force of friction.\n")


"""
bullshit formula for a system of two coupled pulleys
Parameters: mass 1 (smaller mass), radius 1 (corresponding with mass 1), mass 2 (larger mass), radius 2 (corresponding with mass 2),
            moment of inertia of system, acceleration (defaults to g)
"""
def alpha_2pulleys(m1, r1, m2, r2, I, a = g):
    return ( m1*r2*a - m2*r2*a) / ( I + m1*r1**2 + m2*r2**2 )


# KINEMATICS


# gives urls to images of formulas
def kinematics_info():

    print("Regular kinematics:")
    print("https://present5.com/presentation/d3a25a6f8043b1e768d6e03fa3840b6f/image-61.jpg")
    print("Angular kinematics:")
    print("https://i.pinimg.com/736x/35/6f/91/356f91832530bee12258bf8e55b3bfb9--astrophysics-calculus.jpg\n")



# WORK & ENERGY


"""
Basic work formula 
Solves for whatever's <None>
"""
def work_info(W, F, d):

    # calc whatever's unknown
    if W == None: W = F * d
    elif F == None: F = W / d
    elif d == None: d = W / F

    # print info after calcs
    print("Work:", W)
    print("Force (N):", F)
    print("Displacement (m):", d)

    return (W, F, d)


# KE and U functions?
# prob not, conservation of energy gets complicated & it would take too long to implement


# gives basic formulas & explains conservation of energy
def energy_help():

    print("ENERGY")
    print("\tKinetic energy (KE): (1/2)mv^2")
    print("\tPotential energy (U): mah\tmass * acceleration * height")
    print("\tUnits: Joules")
    print("If energy is conserved, then:")
    print("\tinitial potential energy + initial kinetic energy = final potential energy + final kinetic energy\n")



# MOMENTUM


"""
Finds basic momentum info based on what parameter is <None>
Parameters: momentum, mass (kg), velocity (m/s) 
"""
def p_info(p, m, v):
    
    # calc what's not known
    if p == None: p = m * v
    if m == None: m = p / v
    if v == None: v = p / m

    # print results after calcs
    print("Momentum:", p)
    print("Mass (kg):", m)
    print("Velocity (m/s):", v)

    return (p, m, v)


# function to ask oracle for help on momentum formulas
def p_help():

    print("MOMENTUM FORMULAS")
    print("\tp = mv\n")
    print("Elastic collisions:")
    print("\tBoth KE and p are conserved.")
    print("\t(inital p 1) + (initial p 2) = (final p 1) + (initial p 2)")
    print("\tInitial p 1 = final p 2, Initial p 2 = final p 1\n")
    print("Inelastic collisions:")
    print("\tOnly p is conserved.")
    print("(initial p 1) + (initial p 2) = velocity(mass 1 + mass 2)\n")



# CENTRIPEDAL FORCE

"""
Centripedal force info
Parameters: centripedal force (N), mass (kg), velocity (m/s), radius (r)
"""
def f_cp_info(F_cp, m, v, r):

    # whatever's not known, calc it
    if F_cp == None: F_cp = m * v**2 / r
    elif m == None: m = F_cp * r / v**2
    elif v == None: v = sqrt(F_cp * r / m)
    elif r == None: m * v**2 / F_cp

    # centripedal acceleration
    a_cp = v**2 / r

    # print results after calcs
    print("Centripedal force (N):", F_cp)
    print("Centripedal acceleration (m/s**2):", a_cp)
    print("Mass (kg):", m)
    print("Velocity (m/s):", v)
    print("Radius (m):", r)

    return (F_cp, m, v, r)



# calculates the minimum velocity needed to not fall over when spinning in a vertical circle
def min_vertical_v_cp(r):
    return sqrt(g*r)


# function to return info about centripedal force
def f_cp_help():

    print("CENTRIPEDAL FORCE")
    print("Centripedal force formula: F = mv^2 / r")
    print("When spinning in a vertical circle: ")
    print("Bottom of circle: net force = gravity force - normal force")
    print("Top of a circle: net force = gravity force + normal force")




# ANGULAR MOTION & MOMENTUM

# all-encompassing help function for all things angular
def angular_help():

    print("ANGULAR MOTION & MOMENTUM")
    print("The funny Greek letters:")
    print("\talpha (jesus fish) = angular acceleration")
    print("\tomega (w-thing) = angular velocity")
    print("Some helpful conversions, from translational to angular:")
    print("\talpha = a / r")
    print("\tomega = v / r")
    print("Torque formulas:")
    print("\tTorque = (moment of inertia, I)(angular acceleration, alpha)")
    print("\t\tMoment of inertia varies depending on object. For a point particle it's I = mr^2")
    print("\tTorque = (lever arm)(momentum)")
    print("Parallel axis theorem: previous moment of inertia + (mass)(distance to new axis)**2 = new moment of inertia")



# CHAPTER 13: PLANETS AND SHIT

"""
Calculates info about two large (celestial) objects, based on what is left as <None> 
Parameters: force (N), potential energy (J), mass 1 (kg), mass 2 (kg), radius (m)
"""
def grav_force_info(F, U, m1, m2, r):

    # if info relating to F is known, find missing variable & calc U
    if U == None:

        # whatever variable isn't known, calc it
        if F == None: F = G * m1 * m2 / r**2
        elif m1 == None: m1 = (F * r**2) / (G * m2)
        elif m2 == None: m2 = (F * r**2) / (G * m1)
        elif r == None: r = sqrt(G * m1 * m2 / F)
        
        # calc potential energy
        U = G * m1 * m2 / r
    
    # if info related to U is known, find missing variable & calc F
    elif F == None:

        # whatever variable isn't known, calc it
        if U == None: U = G * m1 * m2 / r
        elif m1 == None: m1 = (U * r) / (G * m2)
        elif m2 == None: m2 = (U * r) / (G * m1)
        elif r == None: r = G * m1 * m2 / U

        # calc force
        F = G * m1 * m2 / r**2
    
    # print results after calcs
    print("Force (N):", F)
    print("Potential energy (J):", U)
    print("Mass 1 (kg):", m1)
    print("Mass 2 (kg):", m2)
    print("Radius (m):", r)

    return (F, U, m1, m2, r)


"""
Calculates the missing stat of a planet based on the parameters which aren't <None>
Parameters: acceleration due to gravity on the planet (m/s**2), mass of the planet (kg), radius of the planet (m)
"""
def planet_info(g_planet, M_planet, R_planet):

    # if g isn't given
    if g_planet == None:
        g_planet = G * M_planet / R_planet**2
    
    # if M isn't given
    elif M_planet == None:
        M_planet = g_planet * R_planet**2 / G
    
    # if r isn't given
    elif R_planet == None:
        R_planet = sqrt(G * M_planet / g_planet)
    
    # print results of calcs
    print("Acceleration due to gravity (m/s**2):", g_planet)
    print("Planet mass (kg):", M_planet)
    print("Radius of planet (m):", R_planet)

    return (g_planet, M_planet, R_planet)

"""
Finds escape-velocity related info (like escape velocity itself) based on what parameters aren't <None>
Parameters: escape velocity of the planet (m/s), mass of the planet (kg), radius of the planet (m)
"""
def escape_info(v_escape_planet, M_planet, R_planet):

    # if escape velocity is missing, calc it
    if v_escape_planet == None:
        v_escape_planet = sqrt(2 * G * M_planet / R_planet)
    
    # if planet mass is missing, calc it
    elif M_planet == None:
        M_planet = R_planet * v_escape_planet**2 / (2*G)
    
    # if planet radius is missing, calc it
    elif R_planet == None:
        R_planet = 2*G * M_planet / v_escape_planet**2
    
    # print results of calcs
    print("Escape velocity (m/s):", v_escape_planet)
    print("Planet mass (kg):", M_planet)
    print("Planet radius (m):", R_planet)

    return (v_escape_planet, M_planet, R_planet)


"""
Physics hw 11 problem 7
For determining final velocity of an asteroid from outer space to the atmosphere
Mandatory paremeters: initial velocity (m/s), radius (distance of asteroid from earth, meters)
Optional parameters: radius of the atmosphere (meters), radius of the earth (meters), mass of the earth (kilograms)
"""
def asteroid_v_f(v_i, r, r_atm = 10**5, R_planet = r_earth, M_planet = M_earth):
    return ( v_i**2 - (( 2 * G * M_planet ) / ( r + R_planet ) ) + ((2 * G * M_planet) / (R_planet + r_atm))) ** (1/2)


"""
Physics hw 11 problem 10
Finds 2nd velocity after it's travelled an amount starting at the 1st velocity (v_1, meters)
Mandatory parameters: velocity 1 (m/s), radius 1 (meters), radius 2 (meters)
Optional parameters: mass of the given celestial object (here it's the sun) (kg)
"""
def v_2(v_1, r1, r2, M = 1.98 * 10**30):
    return ( v_1**2 - 2 * G * M * (1/r1 - 1/r2) ) ** (1/2)


"""
HW 11 problem 8 & 9
Finds info of something that's orbiting (involving period)
Parameters: planet period (s), planet mass (kg), planet radius (m)
"""
def planet_period_info(T_planet, M_planet, R_planet):

    # if period is unknown, calc it
    if T_planet == None:
        T_planet = 2*pi * sqrt(R_planet**3 / (G * M_planet))

    # if planet mass is unknown, calc it
    elif M_planet == None:
        M_planet = R_planet**3 / ( G * (T_planet / (2*pi) )**2 )
    
    # if planet radius is unknown, calc it
    elif R_planet == None:
        R_planet = cbrt(G*M_planet * (T_planet/(2*pi))**2)

    # when all said and done
    print("Period (s):", T_planet)
    print("Planet mass:", M_planet)
    print("Planet radius:", R_planet)

    return (T_planet, M_planet, R_planet)



# CHAPTER 14: SIMPLE HARMONIC MOTION

# useful for finding frequency from period, or period from frequency
def reciprocal(n):
    return 1 / n


"""
Finds info about period & related variables in formula, based on what inputs are <None>
Parameters: period (s), mass (kg), k constant
"""
def period_info(T, m, k):

    # if period isn't known, calc it
    if T == None: T = 2*pi * sqrt(m/k)
    
    # if mass isn't known, calc it
    elif m == None: m = k * (T/(2*pi))**2

    # if k const isn't known, calc it
    elif k == None: k = m / (T/(2*pi))**2
    
    # print results of calcs
    print("Period T (s):", T)
    print("Mass m (kg):", m)
    print("k constant:", k)

    return (T, m, k)


"""
Finds info which can be gained from hooke's law, F = -kx
Mandatory parameters: mass (kg), k constant, distance of spring stretched at rest (m)
Optional paramteres: acceleration (m/s**2), defaults to acceleration due to gravity on earth
"""
def hooke_info(m, k, x, a = g):

    # if mass isn't known, calc it
    if m == None: m = -k * x / a

    # if k const isn't known, calc it
    elif k == None: k = m * g / (-x)

    # if x isn't known, calc it
    elif x == None: x = m * g / (-k)

    # if a isn't known (and isn't g), calc it
    elif a == None: a = -k * x / m

    # print results after calcs
    print("Mass (kg):", m)
    print("k constant:", k)
    print("Displacement (m):", x)
    print("Acceleration (m/s**2):", a)
    
    return (m, k, x, a)

"""
Gets info about a spring's potential/kinetic energy, based on what parameters aren't <None>
Parameters: k constant, displacement (m), mass (kg), velocity (m/s), boolean flag to indicate if energy is conserved
"""
def spring_energy_info(k, x, m, v, E_con):

    KE = None
    U = None

    # if all info about potential energy is known, calc it
    if k != None and x != None:

        U = 0.5 * k * x**2

        # if m isn't known but v is, and energy is conserved, solve for m
        if m == None and v != None and E_con: m = U / (0.5 * v**2)

        # if m is't known but v isn't, and energy is conserved, solve for v
        elif v == None and m != None: v = sqrt(U / (0.5 * m))
            
    # if all info about kinetic energy is known and energy is conserved, calc it
    elif m != None and v != None and E_con:

        KE = 0.5 * m * v**2

        # if k isn't known but x is, and energy is conserved, solve for k
        if k == None and x != None: k = KE / (0.5 * x**2)

        # if k is known but x isn't, and energy is conserved, solve for x
        if k != None and x == None: x = sqrt(KE / (0.5 * v))

    # if one type of energy is solved for and energy is conserved, set the other energy equal to it
    if E_con:
        if KE != None and U == None: U = KE
        elif KE == None and U != None: KE = U

    # after all calculations, print all variables
    print("Potential energy (J):", U)
    print("k constant:", k)
    print("Displacement (m):", x)
    print("Kinetic energy (J):", KE)
    print("Mass (kg):", m)
    print("Velocity (m/s):", v)

    return (U, k, x, KE, m, v)


"""
Gets into about a period, but using using variables L (length) and a (acceleration)
instead of m (mass) and k (k constant)
based on what parameters aren't <None>
Mandatory parameters: period (s), length (m)
Optional parameters: a (m/s**2) = acc. due to gravity by default
"""
def period_info_alt(T, L, a = g):

    # if period isn't known, calc it
    if T == None: T = 2*pi * sqrt(L/a)
    # if length isn't known, calc it
    elif L == None: L = a * (T/(2*pi))**2
    # if acc. isn't known (or is default), calc/re-calc it
    elif a == None or a == g: a = L / (T/(2*pi))**2

    # print results after calcs
    print("Period (s):", T)
    print("Length (m):", L)
    print("Acceleration (m/s**2):", a)

    return (T, L, a)



# MAIN

def main():

    # message to confirm module is loaded
    print("\'oracle_help()\' for help")


if __name__ == "main":
    main()
