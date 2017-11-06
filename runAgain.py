import pylab, optparse, os, shutil

#Method that runs the eulerSolve2D program in C It makes a plot of
# the estimated solutions 
#------------------------------------------------------------------
def changeParamEntry(filename, variable,value):
    #Initialize parmeter file
    # variable is the name of the variable in quotes
    # and values is the value of the variable
    #--------------------------------------------------------
    
    pfile = open(filename, 'r')

    
    #Time step
    
    
    lines = pfile.readlines()
    strvals = [l.strip().split() for l in lines]

    inputNameList = list();
    inputValuList = list();
    
    for v in strvals:
        
        inputNameList.append(v[0]);
        if(v[0] == variable):
            inputValuList.append(value)
        else:
            inputValuList.append(v[2]);

    
    
    pfile.close()

    
    pfile = open(filename, 'w')

    for infile_name,infile_value in zip(inputNameList , inputValuList):
        
        infile_value = str(infile_value)
        pfile.write("%s = %s \n" % (infile_name, infile_value) )
        
    pfile.close()

def readParamEntry(filename, variable):
    #Initialize parmeter file
    # variable is the name of the variable in quotes
    # and values is the value of the variable
    #--------------------------------------------------------
    
    pfile = open(filename, 'r')

    
    #Time step
    
    
    lines = pfile.readlines()
    strvals = [l.strip().split() for l in lines]

    inputNameList = list();
    inputValuList = list();
    
    for v in strvals:
        
        inputNameList.append(v[0]);
        if(v[0] == variable):
            pfile.close()
            return v[2]
    pfile.close()
    return 0

def copyResults(dest):
    
    #Danger I am going to remove the existing dest directory
    os.system("rm -r -v "+dest)


    os.system("mkdir "+dest)
    print "copying files..."
    shutil.copy("run/InitialRho.cfg", dest)
    shutil.copy("run/Initialpx.cfg", dest)
    shutil.copy("run/Initialpy.cfg", dest)
    shutil.copy("run/InitialE.cfg", dest)
    shutil.copy("run/Initialpre.cfg", dest)
    shutil.copy("run/Initialvx.cfg", dest)
    shutil.copy("run/Initialvy.cfg", dest)
    shutil.copy("run/Initialx.cfg", dest)
    shutil.copy("run/Initialy.cfg", dest)



    shutil.copy("run/FinalRho.cfg", dest)
    shutil.copy("run/Finalpx.cfg", dest)
    shutil.copy("run/Finalpy.cfg", dest)
    shutil.copy("run/FinalE.cfg", dest)
    shutil.copy("run/Finalpre.cfg", dest)
    shutil.copy("run/Finalvx.cfg", dest)
    shutil.copy("run/Finalvy.cfg", dest)
    shutil.copy("run/Finalx.cfg", dest)
    shutil.copy("run/Finaly.cfg", dest)

    shutil.copy("param.cfg", dest)
    print "done"    

def copyResultsBack(origin):
    
    #Danger I am going to remove the existing origin directory
    #os.system("rm -r -v "+origin)


    #os.system("mkdir "+origin)
    print "copying files..."
    shutil.copy(origin+"/InitialRho.cfg","run/InitialRho.cfg")
    shutil.copy(origin+"/Initialpx.cfg","run/Initialpx.cfg")
    shutil.copy(origin+"/Initialpy.cfg","run/Initialpy.cfg")
    shutil.copy(origin+"/InitialE.cfg","run/InitialE.cfg")
    shutil.copy(origin+"/Initialpre.cfg","run/Initialpre.cfg")
    shutil.copy(origin+"/Initialvx.cfg","run/Initialvx.cfg")
    shutil.copy(origin+"/Initialvy.cfg","run/Initialvy.cfg")
    shutil.copy(origin+"/Initialx.cfg","run/Initialx.cfg")
    shutil.copy(origin+"/Initialy.cfg","run/Initialy.cfg")


    shutil.copy(origin+"/FinalRho.cfg","run/FinalRho.cfg")
    shutil.copy(origin+"/Finalpx.cfg","run/Finalpx.cfg")
    shutil.copy(origin+"/Finalpy.cfg","run/Finalpy.cfg")
    shutil.copy(origin+"/FinalE.cfg","run/FinalE.cfg")
    shutil.copy(origin+"/Finalpre.cfg","run/Finalpre.cfg")
    shutil.copy(origin+"/Finalvx.cfg","run/Finalvx.cfg")
    shutil.copy(origin+"/Finalvy.cfg","run/Finalvy.cfg")
    shutil.copy(origin+"/Finalx.cfg","run/Finalx.cfg")
    shutil.copy(origin+"/Finaly.cfg","run/Finaly.cfg")

    shutil.copy(origin+"/param.cfg","param.cfg")
    print "done"    
        
    
    

#  run c program
#-----------------------------
def runEuler():
    #  run c program with mpi
    #-----------------------------
    #with 2 processos
    #os.system("mpiexec -np 8 ./eulerSolve")
    os.system("./eulerSolve2D")



def plotEulerConsVals(Ufile):#Don't use this
    
    # Plot Results
    #-------------------------------
     


    #Initial Conditions 2nd Order
        
    infile_name = Ufile
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x         = [float(v[0]) for v in strvals]
    y         = [float(v[1]) for v in strvals]
    rho       = [float(v[2]) for v in strvals]
    px        = [float(v[3]) for v in strvals]
    py        = [float(v[4]) for v in strvals]
    E         = [float(v[5]) for v in strvals]

    pylab.figure()
    '''
    pylab.subplot(2,2,1)
    pylab.title("Density")
    pylab.contourf(x,y,rho)
    pylab.title("X Momentum")
    pylab.contourf(x,y,px)
    pylab.title("Y Momentum")
    pylab.contourf(x,y,py)
    pylab.title("Energy")
    pylab.contourf(x,y,E)
    '''
def plotXslice(UfileX,Ufile):
    #Get the x values
    infile_name = UfileX
    infile = open(infile_name)
    lines = infile.readlines()

    x = list()
    strvals = [l.strip().split() for l in lines]
    [x.append(float(v[0])) for v in strvals]

    #Get the array values
     
    infile_name = Ufile
    infile = open(infile_name)
    lines = infile.readlines()

    val = list()
    strvals = [l.strip().split() for l in lines]
    Ny = len(strvals)
    [val.append(float(vi)) for vi in strvals[Ny/2]]
    
    pylab.array(val)
    pylab.plot(x,val,label='X slice')


def plotYslice(UfileY,Ufile):
    #Get the y values
    infile_name = UfileY
    infile = open(infile_name)
    lines = infile.readlines()

    x = list()
    strvals = [l.strip().split() for l in lines]
    [x.append(float(v[0])) for v in strvals]

    #Get the array values
     
    infile_name = Ufile
    infile = open(infile_name)
    lines = infile.readlines()

    val = list()
    strvals = [l.strip().split() for l in lines]
    Ny = len(strvals)
    [val.append(float(v[Ny/2])) for v in strvals]
    
    pylab.array(val)
    pylab.plot(x,val,label='Y slice')

def plot45slice(UfileX,Ufile):
    #Get the y values
    infile_name = UfileX
    infile = open(infile_name)
    lines = infile.readlines()

    x = list()
    strvals = [l.strip().split() for l in lines]
    [x.append(float(v[0])) for v in strvals]
    d = list()
    
    [d.append(pylab.sqrt(2)*xi) for xi in x]

    #Get the array values
     
    infile_name = Ufile
    infile = open(infile_name)
    lines = infile.readlines()

    val = list()
    strvals = [l.strip().split() for l in lines]
    i=0;
    for v in strvals:
        val.append(float(v[i]))
        i = i+1
    
    pylab.plot(d,val,label='Diagonal slice')
    


def plotValue(UfileX, UfileY, Ufile):
    #Get the x values
    infile_name = UfileX
    infile = open(infile_name)
    lines = infile.readlines()

    x = list()
    strvals = [l.strip().split() for l in lines]
    [x.append(float(v[0])) for v in strvals]
    

    #Get the yvalues
    infile_name = UfileY
    infile = open(infile_name)
    lines = infile.readlines()

    y = list()
    strvals = [l.strip().split() for l in lines]
    [y.append(float(v[0])) for v in strvals]
    
    #Get the array values
     
    infile_name = Ufile
    infile = open(infile_name)
    lines = infile.readlines()
    totalMass=0#check to see if we are losing mass
    rho = list()
    strvals = [l.strip().split() for l in lines]
    for v in strvals:
        rhox = list()
        [rhox.append(float(vi)) for vi in v]
        rho.append(rhox)
        for rhoij in rhox:
            totalMass = totalMass+rhoij
        
    print infile_name+" "+str(totalMass)
    
    
    pylab.contourf(x,y,rho)
    pylab.colorbar()

def runEuler1times(dt, title, UfileX, UfileY, Uinitial, Ufinal,
                InitialType):
    changeParamEntry("param.cfg","InitialType",InitialType)
    changeParamEntry("param.cfg","runTime",dt)
    changeParamEntry("time.cfg","t",0)
    time = readParamEntry("time.cfg", "t")
    runEuler()
    
    pylab.figure()
    pylab.suptitle(title)
    pylab.subplot(2,1,1)
    plotValue(UfileX, UfileY, Uinitial)
    pylab.xlabel("x at t=" + time)


    time = readParamEntry("time.cfg", "t")
    pylab.subplot(2,1,2)
    plotValue(UfileX, UfileY, Ufinal)
    pylab.xlabel("x at t= " + time)
    
    
def runEulerNtimes(N, dt, title, InitialType, dest, picDest):
    count = 0
    
    changeParamEntry("param.cfg","InitialType",InitialType)
    changeParamEntry("param.cfg","runTime",dt)
    #this version the time is not reset
    time = readParamEntry("time.cfg", "t")
    strTime = str(time)
    
    runEuler()
    count = count+1
    
    
    numString = str(count)
    #copyResults(dest+"iter"+numString)


    fig1=pylab.figure()
    fig1.subplots_adjust(left=.125, bottom=.1, right=.9, top=.9,
                wspace=.4, hspace=.3)

    pylab.suptitle(title + "at time " + strTime)
    pylab.subplot(2,2,1)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/InitialRho.cfg")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.title("Density")
    pylab.subplot(2,2,2)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Initialvx.cfg")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.title("Px")
    pylab.subplot(2,2,3)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Initialvy.cfg")
    pylab.title("Py")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.subplot(2,2,4)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Initialpre.cfg")
    pylab.title("Pressure")
    pylab.xlabel("x")
    pylab.ylabel("y")
    
    pylab.savefig(picDest+strTime+".png")

    time = readParamEntry("time.cfg", "t")
    strTime = str(time)
    
    fig2=pylab.figure()
    fig2.subplots_adjust(left=.125, bottom=.1, right=.9, top=.9,
                wspace=.4, hspace=.3)



    pylab.suptitle(title + "at time " + strTime)
    pylab.subplot(2,2,1)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/FinalRho.cfg")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.title("Density")
    pylab.subplot(2,2,2)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalvx.cfg")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.title("Px")
    pylab.subplot(2,2,3)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalvy.cfg")
    pylab.title("Py")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.subplot(2,2,4)
    plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalpre.cfg")
    pylab.title("Pressure")
    pylab.xlabel("x")
    pylab.ylabel("y")
    pylab.savefig(picDest+strTime+".png")
   
    
    changeParamEntry("param.cfg","InitialType",0)
    
    for i in range(N):
        runEuler()
        count = count+1
        numString = str(count)
        #copyResults(dest+"iter"+numString)

        time = readParamEntry("time.cfg", "t")
        strTime = str(time)

        fig3=pylab.figure()
        fig3.subplots_adjust(left=.125, bottom=.1, right=.9, top=.9,
                wspace=.4, hspace=.3)

        pylab.suptitle(title + "at time " + strTime)
        
        pylab.subplot(2,2,1)
        plotValue("run/Finalx.cfg","run/Finaly.cfg","run/FinalRho.cfg")
        pylab.xlabel("x")
        pylab.ylabel("y")
        pylab.title("Density")
        pylab.subplot(2,2,2)
        plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalvx.cfg")
        pylab.xlabel("x")
        pylab.ylabel("y")
        pylab.title("Px")
        pylab.subplot(2,2,3)
        plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalvy.cfg")
        pylab.title("Py")
        pylab.xlabel("x")
        pylab.ylabel("y")
        pylab.subplot(2,2,4)
        plotValue("run/Finalx.cfg","run/Finaly.cfg","run/Finalpre.cfg")
        pylab.title("Pressure")
        pylab.xlabel("x")
        pylab.ylabel("y")





        pylab.savefig(picDest+strTime+".png")



#changeParamEntry("param.cfg","InitialType",0)    
#runEuler()






print "ran"
origin = "temp"
copyResultsBack(origin)
    

dt = .3


IniTyp = 0 #0)read in data 1)Uniform, 2)Horizontal Strips, 3)Diagonal Strips,
           #4)Horizontal with
           #pertubation, 5)Inner Square, 6)Inner Circle, 7)Horizonatal Strips
           #pressure for gravity 


Num = 10

runEulerNtimes(Num,dt,"Values ",  IniTyp,"12_17_09dGRho", "results/withgravity/KelvinHemGrav/12_22_09g/")


dest = "temp"
copyResults(dest)
#Kevin Helmholtz instability
#runEuler6times(dt, "Pressure vs t", "run/Finalx.cfg", "run/Finaly.cfg",
  #             "run/FinalPressure.cfg", 3)



#pylab.show()
