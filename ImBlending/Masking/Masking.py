import pylab as pl
import numpy as np
from scipy import misc
from roipoly import roipoly 
from PIL import Image

# ----------------------------------------------------------------
# ----------------------------------------------------------------
# TO MODIFY
# Load images in BW - as enough for masking
img = misc.imread('destination.jpg', 'L' ) # destination
img2 = misc.imread('fg.jpg', 'L' ) # source non rescaled
# ----------------------------------------------------------------
# ----------------------------------------------------------------



# ----------------------------------------------------------------
# ----------------------------------------------------------------
# THIS IS NOT NECESSARY TO MODIFY UNLESS YOU WANT SOMETHINg FANCY

# show the image
pl.imshow(img, interpolation='nearest', cmap="Greys")
pl.title("left click: line segment         right click: close region")

# let user draw  ROI
ROI1 = roipoly(roicolor='r') #let user draw first ROI

# show the image with the  ROI
pl.imshow(img, interpolation='nearest', cmap="Greys")
ROI1.displayROI()

# show the image with ROI and their mean values
pl.imshow(img, interpolation='nearest', cmap="Greys")
ROI1.displayROI() 
pl.title('The mask')
pl.show()

msk=ROI1.getMask(img)
# ----------------------------------------------------------------
# ----------------------------------------------------------------



#positioning

m, n=np.shape(img2) # size of  destination image (assumed larger)
m2, n2=np.shape(img) # size of  source image

# ----------------------------------------------------------------
# ----------------------------------------------------------------
#TO MODIFY


xd=1 # x shift user defined
yd=1 # y shift user defined

# ----------------------------------------------------------------
# ----------------------------------------------------------------

Mask=np.zeros((m ,n))
Mask[yd:m2+yd, xd:n2+xd]=np.array(msk.astype('float'))
pl.imshow(Mask)
pl.show()

# white is 255
Mask=Mask.astype('uint8')*255
im = Image.fromarray(Mask)
im.save("mask.jpg")


# Position source on same domain
fstar = misc.imread('bg.jpg' )*0.
gstar = misc.imread('fg.jpg' )
fstar[yd:m2+yd, xd:n2+xd, :]=gstar
misc.imsave("source.jpg", fstar)









