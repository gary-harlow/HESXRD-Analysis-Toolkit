import mayavi.mlab as mlab
import numpy as np

#this is a temporary solution and should instead be implemneted in qt
def scene(hat):
    
    grid=hat.showData
    mlab.pipeline.volume(mlab.pipeline.scalar_field(grid),vmin=5, vmax=150)
    mlab.outline()
    mlab.show()

