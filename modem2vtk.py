import os
import sys
from modemviz import *

sysarg = sys.argv[1]

# Get input params
with open(sysarg, 'r') as f:
  args = f.readlines()
args = list(map(lambda x: x.split('#')[0].strip(' '), args))
root = os.path.basename(args[0])
path, path_dat, path_mod, stn_flag, m0_flag, interval = args

stn_flag = bool(int(stn_flag))
m0_flag = bool(int(m0_flag)) 

# Setup ouput path
vtkloc = os.path.join(path, 'VTK')
if not os.path.exists(vtkloc):
  os.mkdir(vtkloc)

# Generate VTR resitivity files
pngloc = os.path.join(vtkloc, root)
if interval == 'last':
  it0 = -1; itN = None
elif interval == 'all':
  it0 = 0; itN = None
else:
  try:
    it0, itN = map(int, interval.split(' '))
    itN += 1
  except:
    print("Invalid convert iteration option: try 'all', 'last', '[index0] [indexN]'")
    sys.exit()

generateVTK(path, vtkloc, first_iter=it0, last_iter=itN)  

if stn_flag:
  # Generate station VTP file
  vtkfile = os.path.join(vtkloc, 'stations.vtp')
  if os.path.isfile(path_dat):
    data2vtk(path_dat, vtkfile)

if m0_flag:
  # Generate initial model VTK
  vtk_mod = os.path.join(vtkloc, 'initial')
  model2vtk(path_mod, vtk_mod) 
