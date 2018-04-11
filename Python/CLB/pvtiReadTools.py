#pvtiReadTools.py
from paraview.simple import *
from numpy import zeros

def readSlice(myFile, prms, normal, origin, nx, ny):
    """ 
    ##OVERVIEW:
    Read in a .pvti file with a set of supplied parameters e.g. ['Density','Velocity',...]
    these will be model specific. Also takes in maximum domain of the slice nx, ny as well 
    as the origin and normal direction of the slice.

    ##INPUTS:
    myFile: string consisting of file name .pvti        e.g. "output_VTK_P00_000000.pvti"
    prms  : array of strings, quantities from TCLB      e.g. ['Density','Velocity',...]
    normal: normal to the slice that we want to extract e.g. [0,0,1]
    origin: denotes a location for slice                e.g. [0,0,0]     
    nx    : used to extract slice take from .xml file
    ny    : used to extract slice take from .xml file

    NOTES: 
        - not set up to read in vectors
        - e.g. post script vis: 
            import matplotlib.pyplot as plt
            plt.imshow(data['Variable'])
            plt.show()
    """
    nx +=1
    ny +=1
    f = XMLPartitionedImageDataReader(FileName=myFile, 
                                      CellArrayStatus=prms)
    f = CellDatatoPointData(Input=f)

    mySlice = Slice(Input=f,
                    SliceType='Plane',
                    SliceOffsetValues=[0.0])

    mySlice.SliceType.Origin = origin
    mySlice.SliceType.Normal = normal

    myData = servermanager.Fetch(mySlice)
    myData = myData.GetPointData()

    out = dict()
    for i in range(len(prms)):
        data = myData.GetArray(prms[i])
        tmp = zeros([nx,ny])
        for y in range(ny):
            for x in range(nx):
                tmp[x,y] = data.GetTuple(y*nx+x)[0]

        out[prms[i]] = tmp

    return out
