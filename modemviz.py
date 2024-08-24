import os
import glob
import numpy as np
import pandas as pd
import pyevtk as evtk
import vtk


AIR_VALUE = 39.143

def generateVTK(path, vtkloc, first_iter=0, last_iter=None):
  for i in glob.glob(os.path.join(path, '*_NLCG_*.rho'))[first_iter:last_iter]:
    modfile = i
    modname = os.path.basename(modfile).strip('.rho')
    outfile = os.path.join(vtkloc, modname)
    model2vtk(modfile, outfile)


def model2vtk(modfile, vtkfile):
        
    with open(modfile, 'r') as f:
      modcontent = f.readlines()

    nx, ny, nz, tag0, tag = modcontent[1].strip(' \n').split()
    nx = int(nx); ny = int(ny); nz = int(nz)

    dx = list(map(float, modcontent[2].strip(' \n').split()))
    dy = list(map(float, modcontent[3].strip(' \n').split()))
    dz = list(map(float, modcontent[4].strip(' \n').split()))

    count = 5
    layers = []
    for nl in range(nz):
        count += 1
        layers.append(list(map(lambda x: list(map(float, x.strip(' \n').split())), modcontent[count: count+ny])))
        count = count + ny

    if not modcontent[count].strip(' \n').split():
        count += 1

    x0, y0, z0 = list(map(float, modcontent[count].strip(' \n').split()))
    tag_end = float( modcontent[count+1])

    x = np.ones(nx+1) * x0
    x[1:] = x[1:] + np.cumsum(dx)
    x = -x

    y = np.ones(ny+1) * y0
    y[1:] = y[1:] + np.cumsum(dy)

    z = np.ones(nz+1) * z0
    z[1:] = z[1:] + np.cumsum(dz)

    rho = np.array(layers)
    rho = np.swapaxes(rho, 0, 2)
    rho = np.round(np.exp(rho), 5)
    rho[rho > np.exp(AIR_VALUE)] = np.nan

    evtk.hl.gridToVTK(vtkfile, x, y, z, cellData={"Resistivity": rho})




def read_ModEM_data_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Extract header information
    header = [line[1:].strip() for line in lines if line.startswith(">")]
    sign = -1 if ('-' not in header[1]) else 1

    # Extract data lines (those not starting with "#" or ">")
    data_lines = [line for line in lines if not (line.startswith("#") or line.startswith(">"))]
    
    # Extract data columns and convert to appropriate data types
    data = [line.split() for line in data_lines]
    df = pd.DataFrame(data, columns=["period", "station_name", "latitude", "longitude", "elevation", "northing", "easting", "tensor_component", "real_part", "imaginary_part", "error"])
    df = df.astype({"period": "float64", "latitude": "float64", "longitude": "float64", "elevation": "float64", "northing": "float64", "easting": "float64", "real_part": "float64", "imaginary_part": "float64", "error": "float64"})

    # Create complex impedance and convert units from [mV/km]/[nT] to [V/m]/[T]
    df["Z"] = (df["real_part"] + df["imaginary_part"]*1j)
    if sign < 0:
        df["Z"] = np.conj(df["Z"])

    # Calculate apparent resistivity (ohm meters) and phase (degrees)
    mu = 4 * np.pi * 1e-7
    w = 2 * np.pi / df["period"]

    if header[2] == '[mV/km]/[nT]':
        conv_factor = mu * 1000
    else:
        conv_factor = 1

    df["Z"] = df["Z"] * conv_factor
    df["error"] = df["error"] * conv_factor

    dZ = np.real(np.real(df["error"])/np.sqrt(df["Z"] * np.conj(df["Z"])))
    df["rho"] = 1 / (w * mu) * abs(df["Z"])**2
    df["phase"] = np.angle(df["Z"], deg=True)
    df["rho_error"] = np.abs(2 * df["rho"] * dZ)
    df["phase_error"] = (180 / np.pi) * dZ
    
    return df, header




def data2vtk(datfile, vtkfile):
    # Read datafile
    df, _ = read_ModEM_data_file(datfile)

    # Create a vtkPoints object to hold the 3D coordinates of each station
    vtk_points = vtk.vtkPoints()

    # Create a vtkCellArray to specify the connectivity between points
    cells = vtk.vtkCellArray()

    # Create a dictionary to hold VTK arrays for each column of the dataframe
    vtk_arrays = {}

    for column in df.columns:
        if df[column].dtype == 'object':
            vtk_array = vtk.vtkStringArray()
        else:
            vtk_array = vtk.vtkFloatArray()
        vtk_array.SetName(column)
        vtk_arrays[column] = vtk_array

    # Group the dataframe by station name and generate a VTK point for each unique station
    for station_name, group in df.groupby('station_name'):
        x, y, z = group['northing'].mean(), group['easting'].mean(), group['elevation'].mean()
        point_id = vtk_points.InsertNextPoint(z, x, y)
        cells.InsertNextCell(1)
        cells.InsertCellPoint(point_id)

        for column in group.columns:
            vtk_arrays[column].InsertNextValue(group[column].iloc[0])

    # Create a vtkPolyData object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetVerts(cells)

    # Add the arrays to the polydata object as point data
    for vtk_array in vtk_arrays.values():
        polydata.GetPointData().AddArray(vtk_array)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(vtkfile)
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)
    writer.Write()