import json

import numpy as np
from numpy.linalg import norm
import math
import random

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# All Necessary import modules for the project file to run

"""

Classes:

    CelestialBody  // All objects are physics bodies to be affected by gravity using numerical integration methods
        - Satellite (Sub-Class)  // Specialised objects for experiment 3 with added launch parameters
    
    System // The Class holding the simulation and animation code, contains all bodies in a list

Global Attributes:

    experiment (int) : the current experiment being performed
    G (float) : Gravitational constant from file [m^3/kgs^2)
    celestialbodies (list: objects of Celestialbody): the list of all celestialbodies
    parameters_solar (JSON string): the data from the json file
    
    
    {experiment 3 only}
    
        minvel (float): the minimum launch velocity of a satellite [m/s]
        maxvel (float): the maximum launch velocity of a satellite [m/s]
        velinc (float): the increment launch velocity of a satellite [m/s]
        
        minlaunchangle (float): the minimum launch angle of a satellite [degrees]
        maxlaunchangle (float):the minimum launch angle of a satellite [degrees]
        angleinc (float): the increment launch angle of a satellite [degrees]
    
    {experiment 4 only}
    
        planetstoinclude (int): the number of planets out from the sun to include in potential allignment checks
        anglelimit (float): the maximum angle from the average that a planet can be to be considered alligned [degrees]
        
Non-Class Functions:
    
    readfile(filename) return (parameters_solar['num_iterations'],parameters_solar['ageofsim'],parameters_solar['timestep'])
    // Reads attributes of bodies from .json and calls object creation with said attributes
    
    createsystem(parameters_solar) return (parameters_solar['num_iterations'],parameters_solar['ageofsim'],parameters_solar['timestep'])
    // Creates the system object for the simulation to take place
    
    main()
    //  The main function handling user input, shell menu displays and experiment logic
    
    clear()
    // "Clears" the shell
    
    createsatellites(parameters_solar) {experiment 3 only}
    // Creates objects of satellite with the global variables for launch parameters
    
    


"""
    
    
    
    
class CelestialBody:
    """This is the class for each body in the system, it stores all values read from the file as an attribute.
        This is then put into the System class to be used and animated.
        
        
    
    Attributes:
        
            name (str): Name of the celestial body.
            mass (float): Mass of the body. [inputed in Earth Masses converted into kg]
            orbitalradius (float): Orbital radius in meters. [inputed in AU converted into m]
            solarmass (float): Mass of the solar system's central star. [inputed in Earth Masses converted into kg]
            radius (float): Radius of the body in meters. [m]
            colour (str): Color for visualization.
            actualperiod (float): Actual orbital period in seconds from NASA file. [s]
            lastperiod (float): Used to calculate period from current timeelapsed [s]
            pos (numpy float vector array): the current position of the body in 2d space [m]
            time (float): timestep for the program [s]
            vel (numpy float vector array): current velocity of the body [m/s]
            force (numpy float vector array): current force acting on the body [N]
            newacc (numpy float vector array): next acceleration of the body [m/s^2]
            acc (numpy float vector array): current acceleration of the body [m/s^2]
            prevacc (numpy float vector array): previous acceleration of the body [m/s^2]
            poslist (array: numpy float vector array): list of all positions of the body
        
    
    Methods:
    
            __init(timestep,name,mass,orbitalradius,solarmass,radius,colour,actualperiod)__
            // Initialises a CelestialBody with given parameters from the json.
            
            updatepositionBeeman()
            // Updates the position using Beeman's integration method.
    
            updatepositionEueler()
            // Updates the position using Euler's integration method.
            
            updatevelocityBeeman(new_accel) return self.vel
            // Updates the velocity using Beeman's integration method.
            
            updatevelocityEueler(new_accel) return self.vel
            // Updates the velocity using Euler's integration method.
            
            periodcount(timeelapsed)
            // Calculates the time since the last complete orbit storing it as the body's period.
            
            calcKE() return self.KE
            // Calculates the Total Kinetic energy of the system.

"""
    
    def __init__(self,timestep,name,mass,orbitalradius,solarmass,radius,colour,actualperiod):
        """
            Initialises a CelestialBody with given parameters from the json.
            
            input: name,mass,orbitalradius,solarmass,radius,colour,actualperiod
            

        """
        
        
        
        
        
        
        # initialising every variable from the json file
        
        self.actualperiod = actualperiod
        self.name = name
        self.mass = mass
        self.time = timestep
        self.radius = radius
        self.name = name
        self.colour = colour
        
        self.pos = np.array([0,(1.496* 10**11)*orbitalradius])  # initalised as 0 in the x and the orbital radius in the y, converted from AU to m
        self.poslist = []
        
        
        self.newacc = np.array([float(0),float(0)])
        self.acc = self.newacc
        self.prevacc = self.acc # all acceleration values initialised as 0, 0, before the simulation begins they will be calculated to ensure accuracy 
        
        self.force = np.array([float(0),float(0)])
        
        # {experiment 1 only}
        self.lastperiod = 0
        self.periodlist = []
        self.averageperiod = 0
        
        if name != 'sun':
            
            self.vel = np.array([np.sqrt((G*float(solarmass)* 5.978 * (10**24))/float((1.496* 10**11)*orbitalradius)),float(0)])
            # standard velocity formula using the mass of the sun and the orbital radius converted into SI units
        
        else:
            
            self.vel = np.array([float(0),float(0)])
            # sets the suns velocity to 0 to start, though it is later calculated as it is part of the n-body problem
        
        
        
    
    
    def updatepositionBeeman(self):
        """Updates the position using Beeman's integration method."""
        
        
        
        
        self.pos = self.pos + (self.vel * self.time) + ((1/6)*((4*self.acc) - (self.prevacc))*(self.time**2))
        # the Beeman integration formula for position

    
    def updatepositionEueler(self):
        """Updates the Position using Euler's integration method."""
        
        self.pos = self.pos+( self.vel*self.time)
        # the Euler-Crohmer integration formula for position
        
        
        # !!! note self.pos +=  (self.vel*self.time) doesn't function correctly here in Thonny IDE, this may be due to a difference in += and = in numpy arrays.
    
    def updatevelocityBeeman(self,new_accel):
        """Updates the velocity using Beeman's integration method.
        
        takes new_accel and returns self.vel


        """
        
        self.vel = self.vel + ((1/6)*((2*new_accel)+(5*self.acc)-(self.prevacc)) * self.time)
        
        
        # the Beeman integration formula for velocity
        
        
        return self.vel
      
    def updatevelocityEueler(self,new_accel):
        """
        Updates the velocity using Euler's integration method.

        takes new_accel and returns self.vel

        """
        
        
        self.vel = self.vel +( new_accel*self.time)
        # the Euler-Crohmer integration formula for velocity
        
        
        # !!! note self.vel += (new_accel*self.time) doesn't function correctly here in Thonny IDE, this may be due to a difference in += and = in numpy arrays.
        
        return self.vel
    
    def periodcount(self,timeelapsed):
        """
        Calculates the time since the last complete orbit storing it as the body's period.

        input: timeelapsed


        """
        
        
        self.period = timeelapsed - self.lastperiod
        self.lastperiod = timeelapsed
    
    def calcKE(self):
        """
        Calculates the Total Kinetic energy of the system.

        
        returns: self.KE


        """
        
        self.KE = 0.5*(norm(self.vel)**2)*self.mass
        
        # The standard formula for KE = 1/2mv^2
        
        return self.KE
    
    
    
class Satellite(CelestialBody):
    """
        Creates a Satellite with given parameters from the json. Inherits Properties from parent class CelestialBody
            
        Attributes inherited from Parent-class CelestialBody:
        
            name (str): Name of the celestial body.
            mass (float): Mass of the body. [inputed in Earth Masses converted into kg]
            orbitalradius (float): Orbital radius in meters. [inputed in AU converted into m]
            solarmass (float): Mass of the solar system's central star. [inputed in Earth Masses converted into kg]
            radius (float): Radius of the body in meters. [m]
            colour (str): Color for visualization.
            actualperiod (float): Actual orbital period in seconds from NASA file. [s]
            lastperiod (float): Used to calculate period from current timeelapsed [s]
            pos (numpy float vector array): the current position of the body in 2d space [m]
            time (float): timestep for the program [s]
            vel (numpy float vector array): current velocity of the body [m/s]
            newacc (numpy float vector array): next acceleration of the body [m/s^2]
            acc (numpy float vector array): current acceleration of the body [m/s^2]
            prevacc (numpy float vector array): previous acceleration of the body [m/s^2]
            
        New Attributes for Sub-class Satellite:
        
            launchvelocity (numpy float vector array): initial launch velocity of body [m/s]
            launchangle (float): initial launch angle of body [inputed in degrees converted into radians]
            
            
            
        Methods inherited from Parent class CelestialBody:
            
            updatepositionBeeman()
            // Updates the position using Beeman's integration method.
    
            updatepositionEueler()
            // Updates the position using Euler's integration method.
            
            updatevelocityBeeman(new_accel) return self.vel
            // Updates the velocity using Beeman's integration method.
            
            updatevelocityEueler(new_accel) return self.vel
            // Updates the velocity using Euler's integration method.
            
            periodcount(timeelapsed) {not useful for satellite in current structure}
            // Calculates the time since the last complete orbit storing it as the body's period.
            
            calcKE() return self.KE  {not useful for satellite in current structure}
            // Calculates the Total Kinetic energy of the system.
            

        New Methods unique to Satellite's:
        
            __init(timestep,name,mass,orbitalradius,solarmass,radius,colour,actualperiod,launchvelocity,launchangle)__
            // Initializes a Satellite with given parameters.


        """
    def __init__(self,timestep,name,mass,orbitalradius,solarmass,radius,colour,actualperiod,launchvelocity,launchangle):
        """
            Initializes a Satellite with given parameters.
            
            input: name,mass,orbitalradius,solarmass,radius,colour,actualperiod,launchvelocity,launchangle


        """
        
        
        
        
        # initialising every variable from the json file
        self.approachtime = 0
        self.landed = 0 # these will be overwritten later
        
        
        colourlist = ["red","blue","green","yellow","orange","lightblue","green","white","gray","pink","purple","xkcd:sky blue","xkcd:dark green","xkcd:tan","xkcd:salmon","xkcd:hot pink"] # basic colour's needed for satellite depiction
        self.colour = colourlist[random.randint(0,len(colourlist)-1)] # random assignment of colour on plot, if this feature is removed import random is not needed
        
        self.launchvelocity = launchvelocity
        
        self.launchangle = (launchangle*math.pi/180) # converted from radians to degrees
        
        earthvel = [np.sqrt((G*float(solarmass)* 5.978 * (10**24))/float((1.496* 10**11)*orbitalradius)),0] # standard velocity equation Gm/r converted into SI units
        launchvector = np.array([math.cos(self.launchangle)*launchvelocity,math.sin(self.launchangle)*launchvelocity]) # Calculated using trigonometry where x = vel*cos(theta) and y = vel*sin(theta)
        
        self.pos = np.array([0,((1.496* 10**11)*orbitalradius)]) # initalised as 0 in the x and the orbital radius of Earth in the y + a differentiating factor to avoid a div 0 error, converted from AU to m

        self.vel = launchvector + earthvel # velocity calculated by combining launch with the relative motion of the Earth
        
        
        
        
        # ALL BELOW INHERITED FROM PARENT CLASS
        
        self.name = name
        self.mass = mass
        self.time = timestep
        self.radius = radius
        self.name = name
        
        self.newacc = np.array([float(0),float(0)])
        self.acc = self.newacc
        self.prevacc = self.acc # all acceleration values initialised as 0, 0, before the simulation begins they will be calculated to ensure accuracy 
        
        self.poslist = []
        
        
        
        
        
        
        
            

        
        
        
        
        
        
        
    
    
class System:
    """
        This is the class for The whole system used to simulate and animate,
        it stores all other objects in lists. Aswell as its own attributes
        
        
    
    Attributes:
            
            Simulation:
                bodies (list: Objects of CelestialBody and sub-class Satellite): the array storing the bodies for the system
                time  (float): the timestep for the simulation [s]
                timeelapsed (float): the age of the simulation [s]
                integrationtype (str): the method used for the intergation in sim
                numiterations (int): the number of loops to complete in simulation      
                currentinput (str): input by the user
                r (numpy float vector array): displacement vector between 2 bodies [m]
                rhat (numpy float vector array): unit vector for r
                sumofforces (numpy float vector array) : sum of all force vectors from each body [N]
                
                difference (float): the difference in calculated period and actual [s] {experiment 1 only}
                
                elist (array: floats) : array of e values for system {experiment 2 only}
                timelist (array: float) : array of the ages of sim for e graph {experiment 2 only}
                
                closestapproach (float): the closest distance achieved by a satellite to Mars [m] {experiment 3 only} 
                nearestsatellite (Object of Satellite) : the satellite which achieved the closest approach {experiment 3 only} 

                theta (float): angle between sun and body [degrees] {experiment 4 only} 
                thetalist (list: floats): holds all angle measurements of planets [degrees] {experiment 4 only} 
                thetasum (float): sum of all angle measurements of planets [degrees] {experiment 4 only} 
                thetaaverage (float): average value of angles between bodies [degrees] {experiment 4 only} 
                allignmentcount (int): count of the allignments of the planet's {experiment 4 only}
                alligned (boolean) : whether the planets are currently alligned, used to stop multiple counts of the same allignment {experiment 4 only}
                
            Animation:
                timeecount (float) : the age of the sim at frame of the animation [years]
                starscalefactor (int): multiplied with stars radius to make it more visible in the animation
                planetscalefactor (int): multiplied with planets radius to make it more visible in the animation
                satellitescale (int): the width and height of the satellite in the animation exagerrated for viewing ease [m]
                clickid, keyid, closeid (eventID): used to assign methods to events from the animation                          
                fig (plot figure) : the figure shown in the animation
                ax (plot figure axes) : the axes displayed in the animation
                tracking (int): the number of the current body being tracked
                trackingname (string) : the name of the current body being tracked
                windowlifetime (float) : the number of frames the figure has been shown for
                trackingxoffset (float): the x offset provided by the position of the tracked body
                trackingyoffset (float): the y offset provided by the position of the tracked body
                xoffset (float): the x offset provided by the position from the user input
                yoffset (float): the y offset provided by the position from the user input
                zoomscale (float): the zoom variable affected by the user input
                patches (array: patch of graph plot): the artist objects to be applied to the graph in animation
                scalelimit (float): the limit for the 2d axis in space [m]
                
    Methods:
    
            
    
            __init__(bodies,time)
            
            Simulation:
            
                simulate(numiterations)
            
                tick()
            
                distancecalc(body,otherbody)
            
                calculateforce(body,l)
            
                calculateaccel(body,l) return(body.newacc)
            
            Animation:
            
                displayKEgraph(scheme)
            
                displayanimation()
            
                animate(i)
                
                Run Time Event Handling:
                
                    onclick(event)
                    
                    onkeypress(event)
                    
                    onclose(event)
            
    """
    
    
    def __init__(self,bodies,time):
        
        # initialise the System with objects of Celestialbody and initialise variables
        self.bodies = bodies
        self.time = time
        self.timeelapsed = 0
        self.integrationtype = ""
        
        # {experiment 2 only}
        if experiment == 2:
        
            self.elist=[]
            self.timelist=[] # data for plot
        
        # {experiment 3 only}
        elif experiment == 3:
        
            self.closestapproach = (10**100) # arbitrary large number to be overwritten
        
        # (experiment 4 only}
        elif experiment == 4:
            
            self.allignmentcount = 0
            self.alligned = False # initialiseplanetray allignment conditions
        
        
        # initialise the accelerations for each body ensuring they start with a non-zero acceleration for new current and previous
        for i in range(0,3):
            
            for l, body in enumerate(bodies):
                
                self.calculateaccel(body,l)
                
    
    def simulate(self,numiterations):
        """
        The Physics Simulation takes place in this method. this loops over every body for
        every timestep storing their positions for later display, it calls tick() to step
        the sim forwards one timestep.

        """
        f = open("Output.txt", "w")
        f.write("Energies Output: \n")
        f.close()
        
        self.allignmenttimes = []
        self.numiterations = numiterations
        self.timeelapsed= 0
        self.timelist = []
        self.elist=[]
        
        # run the full simulation storing the positions of the bodies in an array
        
        for loop in range(0,numiterations): # loop every timestep
            
            for body in self.bodies: #loop every body
                
                body.poslist.append(body.pos) # store position as 2d Vector
                
                # {experiment 1 only}
                if len(body.poslist)> 10 and experiment == 1:  # begin period counting after a 10 step delay
                    
                    if body.poslist[-1][0]>(self.bodies[0].poslist[-1][0]) and body.poslist[-2][0] < (self.bodies[0].poslist[-1][0]) : # if the body passes the positive y axis
                        
                        body.periodcount(self.timeelapsed) # find period from current time and last period
                        
                        difference = (body.period - (body.actualperiod*365*24*60*60)) # calculate difference between value and NASA page
                        
                        #SHELL output of Period and difference
                        print(body.name+"'s simulated period measured at ",str(body.period)+"s, found @", str(round(self.timeelapsed,2))+ "s in Sim. \n The difference with the actual NASA value is",difference," s \n as a percentage that's",(difference/(body.actualperiod*365*24*60*60)),"\n")
                        body.periodlist.append(body.period)
                        
        
            
            self.tick()  # tick the simulation running through all methods for numerical integration and force/acceleration calculations
        
        if experiment == 1:
            
            for body in self.bodies:
                
                totalperiod = 0
                
                for period in body.periodlist:
                    
                    totalperiod += int(period)
                    
                if len(body.periodlist)>0:
                    
                    
                    body.averageperiod = totalperiod / len(body.periodlist)
                    difference = (body.averageperiod - (body.actualperiod*365*24*60*60))
                    print("Average Period of ",body.name, "is", body.averageperiod, "s The difference with the actual NASA value is",difference," s \n as a percentage that's",(difference/(body.actualperiod*365*24*60*60)),"\n")
            print("\n"*3)
    
    def tick(self):
        
        # reset KE from last tick
        TotalKe = 0
        TotalGPE = 0
        
        # {experiment 4 only}
        
        if (experiment == 4) and (self.timeelapsed > (self.time*100)):
            
            #reset thetalist, thetasum and allignment from last tick
            thetalist=[]
            thetasum = 0
            allignment = True
            
            for l, body in enumerate(self.bodies): #loop through every body
                
                if body.name != "sun" and l < planetstoinclude+1: #dont include sun or outer planet's in allignment
                    
                    self.distancecalc(self.bodies[0],body) # calculate distance between sun and body
                    
                    #theta = (math.asin(self.r[1]/norm(self.r)))
                    theta = math.atan2(self.r[1], self.r[0]) * (180/(math.pi))  # calculate angle between body and positive x axis anticlockwise
                    thetasum = thetasum + theta # total angle for average later
                    #print(body.name,theta)
                    
                    thetalist.append(theta) # store angle to check conditions later
                    
                    
            thetaaverage = thetasum/planetstoinclude # calc average theta
            
             
            for theta in thetalist:
                
                if abs(theta - thetaaverage) > anglelimit: # if angle is too large, planets are not alligned
                    
                    allignment = False
                    self.alligned = False
                    break # no need to check other planet's
                
            if allignment == True and self.alligned == False: # if planets are alligned now and weren't previously
                
                self.allignmenttimes.append(self.timeelapsed)
                self.allignmentcount +=1
                
                self.alligned = True # stop double counting of the same allignment period
                
        
        
        
        # check intergation method to use
        if self.integrationtype == "Beeman":
            
            for l, body in enumerate(self.bodies): # loop every body
                
                body.updatepositionBeeman() # calc position for every body using Beeman
              
            for l, body in enumerate(self.bodies): #loop every body again
                
                # {experiment 2 only}
                if experiment ==2 :
                    
                    TotalKe = TotalKe + body.calcKE()# Sum Ke values
                    
                    for j, otherbody in enumerate(self.bodies):
                        
                        if l != j:
                            
                            self.distancecalc(body,otherbody)
                            
                            GPE = (G*otherbody.mass*body.mass)/norm(self.r)
                            
                            TotalGPE = TotalGPE + GPE
                    
                    
                body.updatevelocityBeeman(self.calculateaccel(body,l)) # calc acceleration and then velocity for every body using Beeman
            
        
            
        elif self.integrationtype == "Eueler":
            
            
            for l, body in enumerate(self.bodies): # loop every body
                
                # {experiment 2 only}
                if experiment ==2 :
                    
                    TotalKe = TotalKe + body.calcKE()# Sum Ke values
                    
                    for j, otherbody in enumerate(self.bodies):
                        
                        if l != j:
                            
                            self.distancecalc(body,otherbody)
                            
                            GPE = (G*otherbody.mass*body.mass)/norm(self.r)
                            TotalGPE = TotalGPE + GPE
                            
                body.updatevelocityEueler(self.calculateaccel(body,l)) # calc acceleration and then velocity for every body using Euler-Crohmer
            
            
            for l, body in enumerate(self.bodies): # loop every body again
                
                body.updatepositionEueler() # calc position for every body using Euler-Crohmer
                
                
        elif self.integrationtype == "DirectEueler":
            
            for l, body in enumerate(self.bodies): # loop every body
                
                body.updatepositionEueler() # calc position for every body using Euler-Crohmer
            
            for l, body in enumerate(self.bodies): # loop every body again
                
                # {experiment 2 only}
                if experiment ==2 :
                    
                    TotalKe = TotalKe + body.calcKE() # Sum Ke values
                    
                    for j, otherbody in enumerate(self.bodies):
                        
                        if l != j:
                            
                            self.distancecalc(body,otherbody)
                            
                            GPE = (G*otherbody.mass*body.mass)/norm(self.r)
                            TotalGPE = TotalGPE + GPE
                        
                body.updatevelocityEueler(self.calculateaccel(body,l)) # calc acceleration and then velocity for every body using Euler_Crohmer
            
            
        # {experiment 2 only} 
        if experiment == 2 :
            
            TotalE = TotalKe - (0.5)*TotalGPE
            f = open("Output.txt", "a")
            string = str(f"Total Energy at {self.timeelapsed} s = {TotalE} J \n")
            f.write(string)
            f.close()
            self.elist.append(TotalE) # track total ke to be displayed
            self.timelist.append(self.timeelapsed) # track time to be displayed
            
            
        # {experiment 3 only}
        elif experiment == 3:
            
                for l, body in enumerate(self.bodies): # loop every body
                    
                    if body.name == "Satellite" and body.landed != 1: # check for satellites pre return to earth
                        
                        self.distancecalc(body,self.bodies[4]) # calc distance between mars and satellite
                        normr = abs(norm(self.r) - self.bodies[4].radius ) # get the magnitude of this vector
                        
                        if normr < self.closestapproach: # check for new closest satellite
                            
                            self.closestapproach = normr # update new closest value
                            self.nearestsatellite = body # keep tracking of winning satellite
                            body.approachtime = self.timeelapsed # track time of event
                        
                        
                        
                        if self.timeelapsed> self.time*10: # wait before starting to check landing condiitions
                            
                            self.distancecalc(body,self.bodies[3]) # calc distance between earth and satellite
                            
                            if  abs(norm(self.r)) <  (self.bodies[3].radius*1.5): # if satellite returns to earth
                                
                                body.landed = 1 # track which satellites have landed
                                
                        
                
        
            
        self.timeelapsed += self.time # increment age of sim
        
    def distancecalc(self,body,otherbody):
        """
        Calculates distance between two bodies in sim and produces a unit vector for direction
        """

        self.r = otherbody.pos - body.pos # calculate distance vector
        self.rhat = self.r / norm(self.r) # calculate unit vector
        
        
    def calculateforce(self,body,l):
        
        """
        Calculates Force between given body in sim and all others
        """
        
        
        sumofforces = np.array([0,0]) # clear sum vector
        
        for j in range (0,len(self.bodies)): # loop over every body
            
            

            if self.bodies[j].name != "Satellite": # dont interact with satellites as their mass is arbitrary, saves computing time
                
                if j != l: # dont calculate for self, div 0 error
                    
                    self.distancecalc(body,self.bodies[j]) # calculate distance between bodies
                    sumofforces = sumofforces + (((G * (self.bodies[j].mass) * (body.mass) )/ ((norm(self.r))**2)) * self.rhat) # Total Force = SUM(F1,F2...Fj) Fi = GMmi/r^2
        
        body.force = sumofforces # store force to be used for acceleration
        
        
        
        
        return(body)
    
    def calculateaccel(self,body,l):
        
        """
        Calculates acceleration for given body in sim
        """
        self.calculateforce(body,l)
        #calc force then use the given formula,a = F/M
        
        body.prevacc = body.acc
        body.acc = body.newacc # cycle accelerations to be used in integration scheme
        body.newacc = body.force / body.mass
        
        return(body.newacc) # return acceleration to be used
        
    
    def displayEgraph(self,scheme):
        """
        Displays Total energy over time graph for integration scheme
        """
        # font dictionary
        font1 = {'family':'serif','color':'red','size':20} #large
        font2 = {'family':'serif','color':'darkred','size':15} #small
        
        plt.figure(facecolor='black') # set window background to black
        plt.plot(self.timelist,self.elist, scalex=True, scaley=True,color='red',label='KE', marker='o', linestyle='dashed',linewidth=1, markersize=2)
        # plot time against KE with custom ARGS for style
        
        # label axes and give title with fonts from dictionary
        plt.xlabel("Time [s]", fontdict = font2) 
        plt.ylabel("Energy [J]", fontdict = font2)
        
        plt.title(f"Total Energy for {scheme} Integration Scheme", fontdict = font1)

        ax = plt.gca() # get plot axes to allow for more customisation
        
        # set chosen display colours
        
        ax.set_facecolor("xkcd:mahogany") # plot background
        ax.xaxis.label.set_color('red')   # x axis label    
        ax.yaxis.label.set_color('red')   # y axis label  
        ax.spines['bottom'].set_color('red') # x axis line
        ax.spines['left'].set_color('red')  # y axis line  
        ax.tick_params(axis='x', colors='white')  # x axis increment lines  
        ax.tick_params(axis='y', colors='white')  # y axis increment lines  
        
        
        plt.show() # display plot window
        

       
    def displayanimation(self):
        """
        Opens Graph window for simulation animation and creates all artist objects
        """
        #initialise animation variables:
        
        self.windowlifetime = 0
        self.zoomscale = 1
        self.xoffset=0
        self.yoffset=0
        self.tracking =0
        
        # create patch list for circle artists for celestial bodies and rectangle artists for satellites
        self.patches = []
        
        
        
        # prompt user for graph accuracy level
        print("""Choose an Animation Type:
        [-] 1 : More Accurate Scaling of bodies, Less intuitive to view
        
        or
        
        [-] 2 : Less Accurate Scaling of bodies, Ease of Viewing
        
        [ENTER  (1 or 2) TO CONTINUE]""")
        
        
        currentinput = int(input())
        clear() # clear Shell
        
        
        if currentinput == 1:
            
            starscalefactor = 1
            planetscalefactor = 1
            satellitescale = 8**6  # innacurate satellite scale for viewing purposes
            
        elif currentinput == 2:
            
            starscalefactor = 30
            planetscalefactor = 100
            satellitescale = 8**9 # innacurate scale for all bodies, larger bodies are easier to see zoomed out
        
        
        if experiment == 3:
            
            self.scalelimit = (30 ** 8)*0.8 # start more zoomed in for experiment 3 as view focuses on Mars
            
        else:
            
            self.scalelimit = (30 ** 8)*10 # start zoomed out at the solar system
            
        
        # create plot elements
        self.fig = plt.figure(num='Orbital simulation of Mars and Phobos') # name window with num=""
        self.ax = plt.axes()
        
        self.ax.set_facecolor('xkcd:black') # set background colours to black
        self.fig.set_facecolor('xkcd:black')
        
        
        
        
        self.circles=[]
            
        for j, body in enumerate(celestialbodies): # loop every body
            
            
            if body.name == "sun":
                
                # add a circle artist to the plot for the Sun, apparent radius seperate to allow for seperate scale factor
                self.patches.append(plt.Circle((body.poslist[0]), body.radius*starscalefactor, color=body.colour, animated=True))
                
            elif body.name == "Satellite":
                
                # add a rectangle artist for every satellite body, size set by variable (innacurate for viewing purposes)
                self.patches.append(plt.Rectangle((body.poslist[0]), width=satellitescale, height=satellitescale, edgecolor=body.colour, facecolor="white"))
                
            else:
                
                # add a circle artist to the plot for the celestial bodies, apparent size multiplied by scale factor
                self.patches.append(plt.Circle((body.poslist[0]), body.radius*planetscalefactor, color=body.colour, animated=True))
                
        
        
        for i in range(0, len(self.patches)): # loop every artist
            
            self.ax.add_patch(self.patches[i]) # add the artist patches to the plot
        
        
        
        
        # set up the axes with labels and colour
        self.ax.axis("scaled") # allow scaling of axis
        self.ax.set_xlim(-self.scalelimit, self.scalelimit) # set x axis limit  !!! This will be updated in animation
        self.ax.set_ylim(-self.scalelimit, self.scalelimit) # set y axis limit  !!! This will be updated in animation
        self.ax.set_xlabel("x coords [m]") # set x axis label
        self.ax.set_ylabel("y coords [m]") # set y axis label
        self.ax.xaxis.label.set_color('red') # set x axis label colour      
        self.ax.yaxis.label.set_color('red') # set y axis label colour
        self.ax.spines['bottom'].set_color('red') # set x axis division line colours
        self.ax.spines['left'].set_color('red') # set y axis division linecolours   
        self.ax.tick_params(axis='x', colors='white') # set x axis division number colours   
        self.ax.tick_params(axis='y', colors='white') # set y axis division number colours
        
        # output to shell 
        print("Simulation Complete! Displaying Animation:")
        
        
        # create the animator
        self.anim = FuncAnimation(self.fig, self.animate, frames=len(self.bodies[0].poslist), repeat=True, interval=0, blit=True)
        

        # create event handlers
        clickid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        keyid = self.fig.canvas.mpl_connect('key_press_event', self.onkeypress)
        closeid = self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        
        # show the plot
        plt.show()
        
        
    def animate(self, i):
        """
        Updates Graphing window for simulation animation based on user input and frame number, aswell as position of bodies
        """
        
        
        
        
        plt.draw() # redraw plot every frame, issues with blitting in pyplot. reduces efficiency but allows for more control
        
        for j, patch in enumerate(self.patches):  #loop every artist patch
            
            if self.bodies[j].name == "Satellite":
                
                patch.set_xy(self.bodies[j].poslist[i]) # satellites are rectangles so have different args
            
            patch.center = (self.bodies[j].poslist[i]) # update body position based on simulation data
       
        
       
        for text in self.ax.texts: # loop through all text elements
            
            text.set_visible(False) # disable text boxes from previous frame. less efficient as disabling and re-creation takes place every frame, but allows for changing elements such as time elapsed and impact is minimal on test machine
            
        # changing variables to display in text
        timecount = round(((self.time*i)/(60*60*24*365)),2) # time elapsed in years in aniamtion
        trackingname = self.bodies[self.tracking].name # current body being tracked
        
        
        
        if trackingname == "Satellite": # Satellites have extra information displayed with launch conditions
        
            trackingname = trackingname + " Launched at: " + str(self.bodies[self.tracking].launchvelocity) + " m/s " + str((self.bodies[self.tracking].launchangle)*(180/(math.pi))) +" degrees"
            # update tracking name variable with extra information
        
        if experiment == 1:
            
            trackingname = trackingname + " Period Measured at " + str(self.bodies[self.tracking].averageperiod) + " s " 

        
        
        
        
        
        
        self.ax.text(0.5, 0.05, f"Solar System Experiment {experiment}", bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes, ha="center", fontstyle="italic")
        # title for animation, bottom middle of screen
        
        self.ax.text(-0.2, 1.100, f"Camera Tracking {trackingname}", bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes, ha="left")
        # tracking information, top left corner of screen
        
        self.ax.text(1.2, 0.5, f"Time Elapsed {timecount} years", bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes, ha="left")
        # age of sim in animation, right middle of screen

        if self.windowlifetime < 500: # deletes control help popup after a given time
                    
                    self.ax.text(0, -0.1, f"Press I & O keys to zoom in and out. Arrow keys to pan camera", bbox={'facecolor': 'white', 'alpha': (1-(0.002*self.windowlifetime)), 'pad': 5}, transform=self.ax.transAxes, ha="center", fontstyle="italic")
                    # input help popup, bottom left corner of screen
                    
                    self.ax.text(-0.2, 1.05, f"[Press T to change target]", bbox={'facecolor': 'red', 'alpha': (0.5-(0.0005*self.windowlifetime)), 'pad': 5}, transform=self.ax.transAxes, ha="left", fontstyle="italic", fontsize='x-small')
                    # target controls popup, top left corner of screen, underneath tracking information
                    
        #{experiment 3 only}            
        if experiment == 3:
            
            self.ax.text(1.2, 0.3,f"Closest Approach to Mars, \n Velocity: {round(self.nearestsatellite.launchvelocity,4)} m/s, \n Angle: {round((self.nearestsatellite.launchangle*180)/math.pi, 4)} degrees \n at time {round((self.nearestsatellite.approachtime/(60*60*24*365)),2)} years \n at {round(self.closestapproach,2)} m above the surface",bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 5},transform=self.ax.transAxes, ha="left")
            # satellite approach popup, right middle of screen
            
            if self.nearestsatellite.landed == 1:
                
                self.ax.text(1.2, 0.1,f" This satellite returned to Earth",bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 5},transform=self.ax.transAxes, ha="left")
                # satellite return popup, right middle of screen
                
            else:
                
                self.ax.text(1.2, 0.1,f" This Satellite did not return to Earth",bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 5},transform=self.ax.transAxes, ha="left")
                # satellite return popup, right middle of screen
        
        #{experiment 4 only} 
        if experiment == 4:
            
            self.ax.text(1.2, 0.3,f"First {round(planetstoinclude)} planets alligned {round(self.allignmentcount)} times \n within {anglelimit} degrees in {round((self.time*self.numiterations)/(60*60*24*365),2)} years",bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 5},transform=self.ax.transAxes, ha="left")
            # allignment popup, right middle of screen
            
            allignmentoccuring = False
            for time in self.allignmenttimes:
                
                if abs((timecount - (time/(365*60*60*24)))) < 0.05:
                    
                    allignmentoccuring = True
                    
            if allignmentoccuring == True:
                
                self.ax.text(1.2, 0.2, f"Allignment occuring", bbox={'facecolor': 'red', 'alpha': 0.8, 'pad': 5}, transform=self.ax.transAxes, ha="left", fontstyle="italic")
                
                
            else:
                
                self.ax.text(1.2, 0.2, f"Allignment not occuring", bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5}, transform=self.ax.transAxes, ha="left", fontstyle="italic")
                
                
        self.trackingxoffset = self.bodies[self.tracking].poslist[i][0]
        self.trackingyoffset = self.bodies[self.tracking].poslist[i][1]
        # set tracking variables based on bodies position
        
        self.ax.set_xlim(((-self.scalelimit+self.xoffset)*self.zoomscale)+self.trackingxoffset, ((self.scalelimit+self.xoffset)*self.zoomscale)+self.trackingxoffset)
        self.ax.set_ylim(((-self.scalelimit+self.yoffset)*self.zoomscale)+self.trackingyoffset, ((self.scalelimit+self.yoffset)*self.zoomscale)+self.trackingyoffset)
        # adjust sclae of graph to move the camera view
        
        
        
        self.windowlifetime += 1 # increment window lifetime, seperate to frame number as it persists past animation reset
        
        return self.patches # return patches to be animated
        
    def onclick(self,event):
        """
        Event handling method for mouse input
        """
        
        if event.button == 1: # Left click
            
            self.zoomscale -= 0.1 * self.zoomscale # zoom in
            
        elif event.button == 3: # Right click
            
            self.zoomscale += 0.1 * self.zoomscale # zoom out
    
    def onkeypress(self,event):
        """
        Event handling method for keyboard input
        """
        
        if event.key == "i": # i key pressed
            
            self.zoomscale -= 0.1 * self.zoomscale
            
        elif event.key == "o": # o key pressed
            
            self.zoomscale += 0.1 * self.zoomscale
            
        
        elif event.key == "right": # > arrow key pressed
            
            self.xoffset += 0.1*self.scalelimit # pan camera right
            
        elif event.key == "left": # < arrow key pressed
            
            self.xoffset -= 0.1*self.scalelimit # pan camera left
            
        elif event.key == "up": # ^ arrow key pressed
            
            self.yoffset += 0.1*self.scalelimit # pan camera up
            
        elif event.key == "down": # V arrow key pressed
            
            self.yoffset -= 0.1*self.scalelimit # pan camera down
            
        elif event.key == "t": # t key pressed
            
            if self.tracking > len(self.bodies) - 2: # loop back to start of bodies list
                
                self.tracking = 0 # reset which body to track
            
            else:
                
                self.tracking += 1 # increment which body to track
            
            # reset keyboard input panning to focus back on a body
            self.xoffset = 0 
            self.yoffset = 0
    
    def on_close(self, event):
        """
        Event handling method for window closed
        """
        
        print('Animation Complete!') # shell output
        
        
       

def readfile(filename):
    """
        Reads attributes of bodies from .json and calls object creation with said attributes
    """
    
    with open(filename) as file: # opens the json file, using a with statement to only execute code while file is open and automatically close it
        
        global parameters_solar
        global celestialbodies
        global G
        
        parameters_solar = json.load(file) # loads the json file to a JSON string
        
        celestialbodies = [] # initial empty array to hold bodies
        
        G = float(parameters_solar["grav_const"]) # Gravitational constant from file
        
        
        for body in parameters_solar['bodies']: # loop through every body in JSON
            
            # Create a celestial body with given parameters from JSON and add the body to celestialbodies list
            celestialbodies.append(CelestialBody(parameters_solar['timestep'],body['name'],body['mass']*5.978 * (10**24),body['orbital_radius'],parameters_solar['bodies'][0]['mass'],body['radius'],body['colour'],body['period']))
        
        #{experiment 3 only}
        if experiment == 3:
            
            # create varied satellites with different conditions
            createsatellites(parameters_solar) 
        
        return createsystem(parameters_solar)


def createsatellites(parameters_solar):
    """
    Creates objects of satellite with the global variables for launch parameters
    """
    
        
    for launchvelocity in np.arange(minvel,maxvel,velinc): # loop over every velocity value
            
        for launchangle in np.arange(minlaunchangle,maxlaunchangle,angleinc): # loop over every angle value
                
            #create satellite and add it to the body list to give to the simulation
            celestialbodies.append(Satellite(parameters_solar['timestep'],"Satellite",1000,parameters_solar["bodies"][3]['orbital_radius']+0.01,parameters_solar['bodies'][0]['mass'],10,"purple",0.8,launchvelocity,launchangle))
                
                
                
                
def createsystem(parameters_solar):
    """
    Creates the system object for the simulation to take place
    """
    
    global SolarSystem # global as to be edited in the main method
    SolarSystem = System(celestialbodies,parameters_solar['timestep']) # create system object
    
    return (parameters_solar['num_iterations'],parameters_solar['ageofsim'],parameters_solar['timestep'])
        
        
        
def clear():
    """
    'Clears' the shell by printing empty lines
    """
    
    print("\n"*100) # print 100 empty lines



def main():
    """
    The main function handling user input, shell menu displays and experiment logic
    """
    
    #display 1st Menu page 
    print (r"""

                               .    '    	               .  "   '
                      .  .  .            	     '    	  '
          "`       .   .	
                                               '     '
        .    '      _______________
                ==c(___(o(______(_()
                            \=\
                            )=\
                          //|\\
                       // ||  \\
                      // ||    \\
                     //  ||      \\
                    //             \\
                    
                        
     
    SOLAR SYSTEM SIMULATION:
    
    This is a simulation of the planets in our Solar System. All data for planets
    is entered via parameters_solar.json in this file.

    Choose an experiment to model in our Solar System.
    [PRESS ENTER TO CONTINUE]""")
    
    input() # Pause for enter key
    
    clear() # clear shell
    
    
    
    #display 2nd menu page
    print("""
    
    Enter 0 for just the simulation of the system, or;
    
     [-] Experiment 1: Measure orbital period of each planet and compare to NASA values in file
    
     [-] Experiment 2: Measure difference in Energy Conservation using Beeman, Euler - Crohmer or Direct-Euler Intergration
     
     [-] Experiment 3: Launch a Satellite from Earth to Mars at the start of the simulation
     
     [-] Experiment 4: Count number of astronomical allignements in the innermost planets

    [ENTER  (0,1,2,3 or 4) TO CONTINUE]
    """)
    
    global experiment
    experiment = int(input()) #take user input
    
    
    clear() # clear shell
    
    
    
    
    # get extra inputs for experiment 3 only
    if experiment == 3:
        
        global minvel
        global maxvel
        global velinc
        global minlaunchangle
        global maxlaunchangle
        global angleinc
        #global variables to be used by createsatellites()
        
        # display menu
        print("""Would you like to proceed to launch with:

        [-] 1 : Suggested Values for velocity and angle
        
        or
        
        [-] 2 : Custom Values for velocity and angle
        
        [ENTER  (1 or 2) TO CONTINUE]""")
        
        
        currentinput = int(input()) # take input
        
        clear() # clear shell
        
        
        if currentinput == 1: # pre-set values found through manual iterative process
            
            minvel = 4339
            maxvel = 4340
            velinc = 0.1
            maxlaunchangle = 33.31
            minlaunchangle = 33.29
            angleinc = 0.001
            
        else: # custom values entered by user
            
            print("what is the minimum velocity? [m/s]")
            
            minvel = float(input())
            clear()
            
            print("what is the maximum velocity? [m/s]")
            
            maxvel = float(input())
            clear()
            
            print("at what increment? [m/s]")
            
            velinc = float(input())
            clear()
            
            print("what is the minimum launch angle? [degrees, where 0 matches the tangential velocity of the Earth]")
            
            minlaunchangle = float(input())
            clear()
            
            print("what is the maximum launch angle? [degrees, where 0 matches the tangential velocity of the Earth]")
            
            maxlaunchangle = float(input())
            clear()
            
            print("at what increment? [degrees]")
            
            angleinc = float(input())
            clear()
        
        
    
    
    
    # text menu display    
    print("""Should the simulation end when:
    [-] 1 : The max number of iterations has been exceeded
    
    or
    
    [-] 2 : The max age of the simulation has been exceeded
    
    [ENTER  (1 or 2) TO CONTINUE]""")
    
    
    currentinput = input() # take input
    
    clear() # clear shell
    
    numiterations, maxtimeelapsed, timestep = (readfile("parameters_solar.json")) # the file name to be read
    
    numiterations = int(numiterations)
    maxtimeelapsed = float(maxtimeelapsed)
    timestep = float(timestep)
       
    if currentinput == "2": # num of iterations from max time elapsed in file
        
        numiterations = int(maxtimeelapsed / timestep)
         
        
          
        
    if experiment == 0: # basic simulation only
        
        SolarSystem.integrationtype = "Eueler"
        print("CALCULATING POSITION'S:")
        SolarSystem.simulate(numiterations)
        SolarSystem.displayanimation()
        
        
    elif experiment == 1: # period tracking
        
        SolarSystem.integrationtype = "Beeman"
        print("TRACKING PERIODS:")
        SolarSystem.simulate(numiterations)
        SolarSystem.displayanimation()
        
        
    elif experiment == 2: # energy calculating and KE graphs
        
        
        print("CALCULATING ENERGIES FOR BEEMAN INTEGRATION:")
        
        SolarSystem.integrationtype = "Beeman"
        SolarSystem.simulate(numiterations) # simulate the solar system using Beeman Integration
        
        
        print("[Press Enter to view Energy graph for Beeman]")
        
        input() # pause for enter key
        
        SolarSystem.displayEgraph("Beeman") # display energy over time
        
        print("[Press Enter to begin Euler-Crohmer Scheme]")
        
        input() # pause for enter key
        
        clear()
        
        print("CALCULATING ENERGIES FOR EULER INTEGRATION:")
        
        SolarSystem.integrationtype = "Eueler"
        SolarSystem.simulate(numiterations) # simulate the solar system using Euler-Crohmer Integration
        
        print("[Press Enter to view Energy graph]")
        
        input()# pause for enter key
        
        SolarSystem.displayEgraph("Euler-Crohmer") # display  energy over time
        
        print("[Press Enter to begin Direct Euler Scheme]")
        
        input()# pause for enter key
        
        clear()
        
        print("CALCULATING ENERGIES FOR DIRECT EULER INTEGRATION:")
        
        SolarSystem.integrationtype = "DirectEueler"
        SolarSystem.simulate(numiterations) # simulate the solar system using Direct Euler Integration
        
        print("[Press Enter to view Energy graph]")
        
        input()# pause for enter key
        
        SolarSystem.displayEgraph("DirectEuler") # display  energy over time
    
    elif experiment == 3: # satellite simulation
        
        
        print("LAUNCHING SATELLITES IN SIMULATION:")
        
        SolarSystem.integrationtype = "Beeman"
        SolarSystem.simulate(numiterations) # simulate physics
        
        SolarSystem.displayanimation() # start animation
        
    elif experiment == 4: # orbital allignment
        
        global planetstoinclude
        global anglelimit
        
        print("How many planets need to allign for this experiment? [Where Mercury, Venus and Earth would be 3]")
        
        planetstoinclude = float(input()) # user input for number of planets to include
        
        clear()
        
        print("What is the maximum angle for allignment? [Degrees]")
        
        anglelimit = float(input()) # user input for angle limit of allignment
        
        clear()
        
        SolarSystem.integrationtype = "Beeman"
        SolarSystem.simulate(numiterations) # simulate physics
        
        SolarSystem.displayanimation() # start animation
        
    
    
        
        
    



main() # Begin the program in the main function