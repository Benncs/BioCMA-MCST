import birem
import birem.vtk 
import vtk 
from read_results import * 

pathres = "/home/benjamin/Documenti/code/cpp/biomc/results/5M.h5"
results = import_results(pathres)

filename="/home/benjamin/Documenti/code/cpp/BIREM_new/out/sanofi/cma_mesh.vtu"





reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(filename)
reader.Update()
grid = reader.GetOutput()
mb = vtk.vtkMultiBlockDataSet()


division_array = vtk.vtkIntArray()
division_array.SetName("division_array")
division_array.SetNumberOfComponents(1)
division_array.InsertNextValue(5)
division_array.InsertNextValue(10)
division_array.InsertNextValue(10)

grid.GetFieldData().AddArray(division_array)

print(grid.GetFieldData())

# for i in range(n_t):
  

#     scalar_array = birem.vtk.mk_scalar(records_p_dist[i],"particle")
#     g1 = reader.GetOutput()
#     g1.GetCellData().AddArray(scalar_array)
#     mb.SetBlock(i, g1)





    # Write the modified grid to a new VTK XML file
writer = vtk.vtkXMLUnstructuredGridWriter()
extension = writer.GetDefaultFileExtension()
writer.SetFileName(f"./results/test_mb.{extension}")
writer.SetInputData(grid)
writer.SetDataModeToBinary()
writer.Update()
writer.Write()