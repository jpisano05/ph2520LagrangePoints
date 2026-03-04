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

    if mu < 0.01:  # use analytical for small mass ratios
        L1 = (1 - (mu/3)**(1/3) - mu) * R
        L2 = (1 + (mu/3)**(1/3) - mu) * R
        L3 = (-1 - 5/12*mu) * R
    else:  # use fsolve for larger mass ratios
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
    #Starting values with Earth/Sun
    m1 = 1.989e30
    m2 = 5.972e24
    R = 1.496e8

    #Tkinter init
    root = Tk()
    root.title("Lagrange Point Simulator")
    root.geometry("1800x900")
    root.configure(bg="gray20")

    #Main frame
    mainframe = Frame(root, bg="gray20", padx=10, pady=10)
    mainframe.grid(column=0, row=0, sticky=(N,W,E,S))
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)

    #Input Frame
    input_frame = ttk.LabelFrame(mainframe, text="Parameters", padding=10)
    input_frame.grid(column=0, row=0, sticky=N+W, padx=10, pady=10)

    #Mass 1
    mass1 = StringVar(value="1.989e30")
    mass1_row = ttk.Frame(input_frame)
    mass1_row.pack(fill=X, pady=2)
    ttk.Label(mass1_row, text="Mass 1:").pack(side=LEFT)
    mass1_entry = ttk.Entry(mass1_row, width=15, textvariable=mass1)
    mass1_entry.pack(side=LEFT, padx=5)
    massUnit = StringVar(value="kg")
    ttk.Label(mass1_row, textvariable=massUnit).pack(side=LEFT)

    #Mass 2
    mass2 = StringVar(value="5.972e24")
    mass2_row = ttk.Frame(input_frame)
    mass2_row.pack(fill=X, pady=2)
    ttk.Label(mass2_row, text="Mass 2:").pack(side=LEFT)
    mass2_entry = ttk.Entry(mass2_row, width=15, textvariable=mass2)
    mass2_entry.pack(side=LEFT, padx=5)
    ttk.Label(mass2_row, textvariable=massUnit).pack(side=LEFT)

    #Distance
    distance = StringVar(value="1.496e8")
    distance_row = ttk.Frame(input_frame)
    distance_row.pack(fill=X, pady=2)
    ttk.Label(distance_row, text="Distance:").pack(side=LEFT)
    distance_entry = ttk.Entry(distance_row, width=15, textvariable=distance)
    distance_entry.pack(side=LEFT, padx=5)
    distanceUnit = StringVar(value="km")
    ttk.Label(distance_row, textvariable=distanceUnit).pack(side=LEFT)

    #Unit Selection
    unit_frame = ttk.Frame(input_frame)
    unit_frame.pack(fill=X, pady=5)
    ttk.Label(unit_frame, text="Mass Unit:").pack(side=LEFT)
    massUnitCombo = ttk.Combobox(unit_frame, values=["kg", "Solar Masses"], state="readonly", textvariable=massUnit, width=12)
    massUnitCombo.pack(side=LEFT, padx=5)

    ttk.Label(unit_frame, text="Distance Unit:").pack(side=LEFT, padx=(10, 0))
    distanceUnitCombo = ttk.Combobox(unit_frame, values=["km", "AU"], state="readonly", textvariable=distanceUnit, width=12)
    distanceUnitCombo.pack(side=LEFT, padx=5)

    #Calculate Button
    ttk.Button(input_frame, text="Calculate",
               command=lambda: calculatePoints(
                   float(mass1.get()) * (1.989e30 if massUnit.get() == "Solar Masses" else 1),
                   float(mass2.get()) * (1.989e30 if massUnit.get() == "Solar Masses" else 1),
                   float(distance.get()) * (1.496e8 if distanceUnit.get() == "AU" else 1),
               )).pack(pady=10)

    #Canvas
    windowDimX = 1200
    windowDimY = 800
    canvas_frame = ttk.Frame(mainframe)
    canvas_frame.grid(column=1, row=0, sticky=N+S+E+W, padx=10, pady=10)
    mainframe.columnconfigure(1, weight=1)
    mainframe.rowconfigure(0, weight=1)

    canvas = Canvas(canvas_frame, width=windowDimX, height=windowDimY, bg="black")
    canvas.pack(fill=BOTH, expand=True)

    #Styling
    style = ttk.Style()
    style.theme_use('clam')
    style.configure("TLabel", foreground="white", background="gray20", font=("Helvetica", 12))
    style.configure("TEntry", font=("Helvetica", 12))
    style.configure("TButton", font=("Helvetica", 12), padding=6)
    input_frame = ttk.LabelFrame(mainframe, text="Parameters", padding=10)
    #input_frame.configure(style="TLabelFrame")

    mass1_entry.focus()
    calculatePoints(m1, m2, R)

    root.mainloop()