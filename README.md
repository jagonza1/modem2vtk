# modem2vtk
Python script to convert ModEM a 3D resistivity model to VTK format for 3D visualization. 

## Convert ModEM output to VTK format

- Create virtual environment from `requirements.txt` file:
  `$ python -m venv myenv`
  `$ myenv\Scripts\activate` (Windows)
  `$ pip install -r requirements.txt`

- Prepare input parameters in `settings.cfg`
  - Path to ModEM inversion files
  - Path to observed data file
  - Path to initial model file
  - Convert stations: '1' True, '0' False 
  - Convert initial model: '1' True, '0' False
  - Convert iterations: 'all', 'last', '0 10' from iteration 0 to iteration 10

- Run `modem2vtk` script:
  `(myenv)$ python modem2vtk settings.cfg`
    
- Check output: written to ModEM inversion filepath, subdirectory VTK
  - [modem_NLCG_iter].vtr : resistivity model
  - initial.vtr : initial model
  - stations.vtp : location of MT sites

- Open and inspect *.vtX files in ParaView
