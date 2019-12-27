#rutherfordExperiment v1.0
#author: Jason Yi
#jtyi@ucdavis.edu
#Description: This program simulates the Rutherford gold goil experiment in 1D
#Outputs position,

from formulaBox import FormulaBox
from pointCharge import PointCharge
import matplotlib.pyplot as plt
import numpy as np

def rutherford():
    #initial start
    fb = FormulaBox()
    
    #record data
    timeTable = []
    posTable = []
    velTable = []
    accTable = []
    keTable = []
    espTable = []
    timeBreak = []
    cumETable= []

    #initial states
    alphaPart = PointCharge((fb.amu * 4),1,(fb.fundamentalCharge * 2),1)
    auNucleus = PointCharge((fb.amu * 197),0,(fb.fundamentalCharge * 79),1)

    auNucleus.setPos1D(10)

    ke = fb.kinetic(alphaPart.getMass(),alphaPart.getVelocity())
    esp = fb.electroPoten(alphaPart.getCharge(),auNucleus.getCharge(),fb.oneDimDist(alphaPart.getPos()[0],auNucleus.getPos()[0]))

    totalEnergy = ke + esp #total energy must be constant
    
    #create true time scale
    for i in range(0,5):
        #create hundreth second time scale, use time this for calculations
        timeAdd = i
        for j in range(0,1000):
            timeAdd += 0.001
            timeBreak.append(timeAdd)
                                    
    #begin simulation. Time start
    for k in range(0,len(timeBreak)):
        #find data from current timestep
        currentTime = timeBreak[k]
        currentPos = alphaPart.getPos()[0]
        currentVel = alphaPart.getVelocity()

        #find current acceleration
        currentAcc = -fb.electroStatic(alphaPart.getCharge(),auNucleus.getCharge(),fb.oneDimDist(alphaPart.getPos()[0],auNucleus.getPos()[0])) / alphaPart.getMass()

        #energy calculations
        currentKE = fb.kinetic(alphaPart.getMass(),currentVel)
        currentESP = fb.electroPoten(alphaPart.getCharge(),auNucleus.getCharge(),fb.oneDimDist(alphaPart.getPos()[0],auNucleus.getPos()[0]))
        currentCumE = currentKE + currentESP
        
        #record data
        timeTable.append(currentTime)
        posTable.append(currentPos)
        velTable.append(currentVel)
        accTable.append(currentAcc)
        keTable.append(currentKE)
        espTable.append(currentESP)
        cumETable.append(currentCumE)

        #timeStep update; acc -> vel -> pos -> ke -> esp -> esf
        #acc does not change throughout timeStep
        timeStepVel = fb.kinomaticTwovf(currentVel,currentAcc,0.01)
        timeStepPos = fb.kinomaticOnexf(currentPos,currentVel,currentAcc,0.01)
        alphaPart.setVelocity(timeStepVel)
        alphaPart.setPos1D(timeStepPos)

        timeStepKE = fb.kinetic(alphaPart.getMass(),timeStepVel)
        timeStepESP = fb.electroPoten(alphaPart.getCharge(),auNucleus.getCharge(),fb.oneDimDist(alphaPart.getPos()[0],auNucleus.getPos()[0]))
        timeStepCumE = timeStepKE + timeStepESP

    #represent data
    plt.figure("Position")
    plt.plot(timeTable,posTable)
    plt.figure("Velocity")
    plt.plot(timeTable,velTable)
    plt.figure("Acceleration")
    plt.plot(timeTable,accTable)
    plt.figure("Energy Diagrams")
    plt.plot(timeTable,keTable,'r')
    plt.plot(timeTable,espTable,'b')
    plt.plot(timeTable,cumETable,'g')
    plt.show()
    return
rutherford()
