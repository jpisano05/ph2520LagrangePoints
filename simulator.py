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
    
    r1 = mu * R
    r2 = (1 - mu) * R

    L1 = fsolve(collinear_lagrange, 0.5 - mu, mu)[0] * R
    L2 = fsolve(collinear_lagrange, 1.1 - mu, mu)[0] * R
    L3 = fsolve(collinear_lagrange, -1.0 - mu, mu)[0] * R

    L4 = [(0.5 - mu) * R, (np.sqrt(3) / 2 * R)]
    L5 = [(0.5 - mu) * R, -(np.sqrt(3) / 2) * R]
    
    print("L1: ", L1)
    print("L2: ", L2)
    print("L3: ", L3)
    print("L4: ", L4)
    print("L5: ", L5)
    print("Body 1 Distance from Barycenter ", r1)
    print("Body 2 Distance from Barycenter ", r2)

    drawTwoBodySystem(m1, m2, R, r1, r2, L1, L2, L3, L4, L5)
    
    return L1, L2, L3

#draws the two body system on the canvas. exact locations aren't to scale, but rather are for visualization purposes
#utilizes scaling to still preserve some aspect of difference in dimensionality
def drawTwoBodySystem(m1, m2, R, r1, r2, L1, L2, L3, L4, L5):
    
    canvas.delete("all")
    canvas.grid(column=4, row=1, rowspan=10, pady=20, sticky=N)
    
    #scaling factors
    #how much greater m1 is than m2
    sizeDiff = max(np.log10(m1/m2), 0.1)
    print(sizeDiff)
    relativeDistance = np.log10(R)
    maxScale = 1
    
    visualSeperation = 400
    scale = visualSeperation / R
    
    #solar mass relative
    defaultRadius = 80
    
    #calculate barycenter position, draw last so its on top
    barycenterX = windowDimX/2
    barycenterY = windowDimY/2
    
    #body 1
    body1Radius = defaultRadius * maxScale
    body1X = barycenterX  - r1 * scale
    body1Y = barycenterY
    
    canvas.create_oval(body1X - (body1Radius/2),
                       body1Y - (body1Radius/2), 
                       body1X + (body1Radius/2),
                       body1Y + (body1Radius/2), 
                       outline="orange", fill="yellow", width=2)
    
    #orbits
    canvas.create_oval(
        barycenterX - r1 * scale,
        barycenterY - r1 * scale,
        barycenterX + r1 * scale,
        barycenterY + r1 * scale,
        outline="white"
    )

    canvas.create_oval(
        barycenterX - r2 * scale,
        barycenterY - r2 * scale,
        barycenterX + r2 * scale,
        barycenterY + r2 * scale,
        outline="white"
    )
    
    #body 2
    body2Radius = defaultRadius / sizeDiff
    if(body2Radius > body1Radius):
        body2Radius = body1Radius
    body2X = barycenterX + r2 * scale
    body2Y = body1Y
    canvas.create_oval(body2X - (body2Radius/2),
                       body2Y - (body2Radius/2), 
                       body2X + (body2Radius/2),
                       body2Y + (body2Radius/2),  
                       outline="green", fill="blue", width=2)
    
    #Lagrange Points
    PointRadius = 10
    L1Position = barycenterX + L1 * scale
    L1Y = body2Y
    canvas.create_oval(L1Position - (PointRadius/2),
                       L1Y - (PointRadius/2), 
                       L1Position + (PointRadius/2),
                       L1Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L2Position = barycenterX + L2 * scale
    L2Y = body2Y
    canvas.create_oval(L2Position - (PointRadius/2),
                       L2Y - (PointRadius/2), 
                       L2Position + (PointRadius/2),
                       L2Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L3Position = barycenterX + L3 * scale
    L3Y = body2Y
    canvas.create_oval(L3Position - (PointRadius/2),
                       L3Y - (PointRadius/2), 
                       L3Position + (PointRadius/2),
                       L3Y + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L4XPosition = barycenterX + L4[0] * scale
    L4YPosition = barycenterY - L4[1] * scale
    canvas.create_oval(L4XPosition - (PointRadius/2),
                       L4YPosition - (PointRadius/2), 
                       L4XPosition + (PointRadius/2),
                       L4YPosition + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    L5XPosition = barycenterX + L5[0] * scale
    L5YPosition = barycenterY - L5[1] * scale
    canvas.create_oval(L5XPosition - (PointRadius/2),
                       L5YPosition - (PointRadius/2), 
                       L5XPosition + (PointRadius/2),
                       L5YPosition + (PointRadius/2),  
                       outline="green", fill="red", width=2)
    
    #draw lines
    #L1 line
    canvas.create_line(L1Position, body2Y, body2X, body2Y, fill="white")
    canvas.create_text(L1Position + (body2X - L1Position) / 2, body2Y + 10, fill="white", text=("{:.3e}".format(abs(L1 - r2) / 1000)) + " km")
    #L2 Line
    canvas.create_line(L2Position, body2Y, body2X, body2Y, fill="white")
    canvas.create_text(L2Position + (body2X - L2Position) / 2, body2Y - 10, fill="white", text=("{:.3e}".format(abs(L2 - r2) / 1000)) + " km")
    #L3 Line
    canvas.create_line(L3Position, body1Y, body1X, body1Y, fill="white")
    canvas.create_text(L3Position + (body1X - L3Position) / 2, body2Y - 10, fill="white", text=("{:.3e}".format(abs(L3 - r1) / 1000)) + " km")
    #L4 Line
    canvas.create_line(L4XPosition, L4YPosition, body2X, body2Y, fill="white")
    canvas.create_line(L4XPosition, L4YPosition, body1X, body1Y, fill="white")
    #L5 Line
    canvas.create_line(L5XPosition, L5YPosition, body2X, body2Y, fill="white")
    canvas.create_line(L5XPosition, L5YPosition, body1X, body1Y, fill="white")
    
    #Draw Barycenter
    canvas.create_line(barycenterX, barycenterY + 10, barycenterX, barycenterY - 10, fill="red")
    canvas.create_line(barycenterX + 10, barycenterY, barycenterX - 10, barycenterY, fill="red")

if __name__ == "__main__":
    #test for earth/sun system
    m1 = 1.989e30
    m2 = 5.972e24
    R = 1.496e8
    
    #pygame init
    root = Tk()
    root.title("Lagrange Point Simulator")
    root.geometry("1800x900")
    
    mainframe = ttk.Frame(root, padding=(3, 3, 12, 12))
    mainframe.grid(column=0, row=0, stick=(N,W,E,S))
    
    windowDimX = 1200
    windowDimY = 800
    canvas = Canvas(mainframe, width=windowDimX, height=windowDimY, bg="black")
    
    mass1 = StringVar(value="1.989e30")
    mass1_entry = ttk.Entry(mainframe, width = 15, textvariable = mass1)
    mass1_entry.grid(column = 2, row=1, sticky=W)
    mass2 = StringVar(value="5.972e24")
    mass2_entry = ttk.Entry(mainframe, width = 15, textvariable = mass2)
    mass2_entry.grid(column = 2, row=2, sticky=W)
    distance = StringVar(value="1.496e8")
    distance_entry = ttk.Entry(mainframe, width = 15, textvariable = distance)
    distance_entry.grid(column = 2, row=3, sticky=W)
    
    ttk.Label(mainframe, text="Mass 1").grid(column=1, row=1, sticky=E)
    ttk.Label(mainframe, text="Mass 2").grid(column=1, row=2, sticky=E)
    ttk.Label(mainframe, text="Distance").grid(column=1, row=3, sticky=E)
    
    ttk.Label(mainframe, text="Unit Selection:").grid(column=1, row=5, sticky=E)
    
    massUnit = StringVar(value="kg")
    distanceUnit = StringVar(value="km")

    massUnitCombo = ttk.Combobox(mainframe, values=["kg", "Solar Masses"], state="readonly", textvariable=massUnit)
    massUnitCombo.grid(column=1, row=6, sticky=E)

    distanceUnitCombo = ttk.Combobox(mainframe, values=["km", "AU"], state="readonly", textvariable=distanceUnit)
    distanceUnitCombo.grid(column=2, row=6, sticky=E)

    massUnitLabel1 = ttk.Label(mainframe, textvariable=massUnit)
    massUnitLabel1.grid(column=3, row=1, sticky=E)

    massUnitLabel2 = ttk.Label(mainframe, textvariable=massUnit)
    massUnitLabel2.grid(column=3, row=2, sticky=E)

    distanceUnitLabel = ttk.Label(mainframe, textvariable=distanceUnit)
    distanceUnitLabel.grid(column=3, row=3, sticky=E)
    
    calculatePoints(m1, m2, R)
    
    ttk.Button(mainframe, text="Calculate",
            command=lambda: calculatePoints(
                float(mass1.get()) * (1.989e30 if massUnit.get() == "Solar Masses" else 1),
                float(mass2.get()) * (1.989e30 if massUnit.get() == "Solar Masses" else 1),
                float(distance.get()) * (1.496e8 if distanceUnit.get() == "AU" else 1),
            )
            ).grid(column=2, row=4, sticky=E)
    
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