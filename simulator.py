import numpy as np
from tkinter import *
from tkinter import ttk
from scipy.optimize import fsolve
import matplotlib

m1 = 0
m2 = 0
R = 0

#helper function for getting the normalized coordinates
def collinear_lagrange(x, mu):
    return x - (1 - mu) * (x + mu)/ abs(x + mu) ** 3 - mu * (x - 1 + mu) / abs(x - 1 + mu) ** 3

#calculate lagrange points then update the drawing
#this code was used as example for this: https://orbital-mechanics.space/the-n-body-problem/Lagrange-points-example.html
def calculatePoints(m1, m2, R):
    R = R * 1000
    
    mu = m2 / (m1 + m2)

    L1 = fsolve(collinear_lagrange, 0.5, mu)[0] * R
    L2 = fsolve(collinear_lagrange, 1.1, mu)[0] * R
    L3 = fsolve(collinear_lagrange, -1.0, mu)[0] * R

    L4 = [(0.5 - mu) * R, (np.sqrt(3) / 2 * R)]
    L5 = [(0.5 - mu) * R, -(np.sqrt(3) / 2) * R]
    
    print("L1: ", L1)
    print("L2: ", L2)
    print("L3: ", L3)
    print("L4: ", L4)
    print("L5: ", L5)


    drawTwoBodySystem(m1, m2, R, L1, L2, L3, L4, L5)
    
    return L1, L2, L3

#draws the two body system on the canvas. exact locations aren't to scale, but rather are for visualization purposes
#utilizes scaling to still preserve some aspect of difference in dimensionality
def drawTwoBodySystem(m1, m2, R, L1, L2, L3, L4, L5):
    windowDimX = 1200
    windowDimY = 800
    canvas = Canvas(mainframe, width=windowDimX, height=windowDimY, bg="black")
    canvas.grid(column=4, row=1, rowspan=10, pady=20, sticky=N)
    
    #scaling factors
    #how much greater m1 is than m2
    sizeDiff = np.log10(m1/m2)
    print(sizeDiff)
    relativeDistance = np.log10(R)
    maxScale = 1

    #solar mass relative
    defaultRadius = 80
    
    #body 1
    body1Radius = defaultRadius * maxScale
    body1X = windowDimX/2
    body1Y = windowDimY/2
    
    canvas.create_oval(body1X - (body1Radius/2),
                       body1Y - (body1Radius/2), 
                       body1X + (body1Radius/2),
                       body1Y + (body1Radius/2), 
                       outline="orange", fill="yellow", width=2)
    
    #orbit
    orbitRadius = (50 * relativeDistance) * maxScale
    canvas.create_oval(body1X - (orbitRadius/2),
                       body1Y - (orbitRadius/2), 
                       body1X + (orbitRadius/2),
                       body1Y + (orbitRadius/2),  
                       outline="white", width=2)
    
    #body 2
    body2Radius = defaultRadius / sizeDiff
    if(body2Radius > body1Radius):
        body2Radius = body1Radius
    body2X = body1X + (orbitRadius / 2)
    body2Y = body1Y
    canvas.create_oval(body2X - (body2Radius/2),
                       body2Y - (body2Radius/2), 
                       body2X + (body2Radius/2),
                       body2Y + (body2Radius/2),  
                       outline="green", fill="blue", width=2)
    
    #Lagrange Points
    PointRadius = 10
    L1Scaled = (L1 / R) * (orbitRadius / 3)
    print(L1Scaled)
    L1Position = body1X + L1Scaled
    L1Y = body2Y
    canvas.create_oval(L1Position - (PointRadius/2),
                       L1Y - (PointRadius/2), 
                       L1Position + (PointRadius/2),
                       L1Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L2Scaled = (L2 / R) * (orbitRadius / 3)
    print(L2Scaled)
    L2Position = body2X + (L2Scaled / 2)
    L2Y = body2Y
    canvas.create_oval(L2Position - (PointRadius/2),
                       L2Y - (PointRadius/2), 
                       L2Position + (PointRadius/2),
                       L2Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L3Scaled = (L3 / R) * (orbitRadius / 2)
    print("L3:")
    print(L3Scaled)
    L3Position = body1X + L3Scaled
    L3Y = body2Y
    canvas.create_oval(L3Position - (PointRadius/2),
                       L3Y - (PointRadius/2), 
                       L3Position + (PointRadius/2),
                       L3Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L4XScaled = (L4[0] / R) * (orbitRadius / 2)
    L4XPosition = L4XScaled + body1X
    L4YScaled = (L4[1] / R) * (orbitRadius / 2)
    L4YPosition = L4YScaled + body1Y
    canvas.create_oval(L4XPosition - (PointRadius/2),
                       L4YPosition - (PointRadius/2), 
                       L4XPosition + (PointRadius/2),
                       L4YPosition + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    print(L4XScaled)
    
    L5XScaled = (L5[0] / R) * (orbitRadius / 2)
    L5XPosition = L5XScaled + body1X
    L5YScaled = (L5[1] / R) * (orbitRadius / 2)
    L5YPosition = L5YScaled + body1Y
    canvas.create_oval(L5XPosition - (PointRadius/2),
                       L5YPosition - (PointRadius/2), 
                       L5XPosition + (PointRadius/2),
                       L5YPosition + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    print(L5XScaled)
    
    #draw lines
   # canvas.create_line()
   
def varUpdate(var1, var2):
    print("Tyrig")
    var1.set(var2)

if __name__ == "__main__":
    #test for earth/sun system
    m1 = 1.989e30
    m2 = 5.972e24
    R = 1.496e11
    
    #pygame init
    root = Tk()
    root.title("Lagrange Point Simulator")
    root.geometry("1200x900")
    
    mainframe = ttk.Frame(root, padding=(3, 3, 12, 12))
    mainframe.grid(column=0, row=0, stick=(N,W,E,S))
    
    calculatePoints(m1, m2, R)
    
    mass1 = StringVar(value="1.989e30")
    mass1_entry = ttk.Entry(mainframe, width = 15, textvariable = mass1)
    mass1_entry.grid(column = 2, row=1, sticky=W)
    mass2 = StringVar(value="5.972e24")
    mass2_entry = ttk.Entry(mainframe, width = 15, textvariable = mass2)
    mass2_entry.grid(column = 2, row=2, sticky=W)
    distance = StringVar(value="1.496e8")
    distance_entry = ttk.Entry(mainframe, width = 15, textvariable = distance)
    distance_entry.grid(column = 2, row=3, sticky=W)
    
    ttk.Button(mainframe, text="Calculate",
            command=lambda: calculatePoints(float(mass1.get()), float(mass2.get()), float(distance.get()))
            ).grid(column=2, row=4, sticky=E)
    
    ttk.Label(mainframe, text="Mass 1").grid(column=1, row=1, sticky=E)
    ttk.Label(mainframe, text="Mass 2").grid(column=1, row=2, sticky=E)
    ttk.Label(mainframe, text="Distance").grid(column=1, row=3, sticky=E)
    
    ttk.Label(mainframe, text="Unit Selection:").grid(column=1, row=5, sticky=E)
    
    massUnits = StringVar(value=["kg", "Solar Masses"])
    massUnitList = Listbox(mainframe, listvariable=massUnits, height=2)
    massUnitList.grid(column=1, row=6, sticky=E)

    distanceUnits = StringVar(value=["km", "AU"])
    distanceUnitList = Listbox(mainframe, listvariable=distanceUnits, height=2)
    distanceUnitList.grid(column=2, row=6, stick=E)
    
    
    massUnit = StringVar(value="kg")
    distanceUnit = StringVar(value="km")
    ttk.Label(mainframe, text=massUnit.get()).grid(column=3, row=1, sticky=E)
    ttk.Label(mainframe, text=massUnit.get()).grid(column=3, row=2, sticky=E)
    ttk.Label(mainframe, text=distanceUnit.get()).grid(column=3, row=3, sticky=E)
    
    print(massUnitList.curselection())
    #massUnit.set(massUnits.get()[massUnitList.curselection()[0]])
    massUnitList.bind('<<ListboxSelect>>', lambda e: print(massUnits.get()))
    
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)	
    #mainframe.columnconfigure(2, weight=1)
    for child in mainframe.winfo_children(): 
        child.grid_configure(padx=5, pady=5)
    mass1_entry.focus()
    
    #style
    style = ttk.Style()
    style.configure("My.TFrame", background="gray")  
    mainframe.configure(style="My.TFrame")
    
    root.mainloop()
    
    #test for earth/sun system
    #m1 = 1.989e30
    #m2 = 5.972e24
    #R = 1.496e11
    
    #L1, L2, L3 = calculatePoints(m1, m2, R)
    
    #print(L1)
    #print(L2)
    #print(L3)