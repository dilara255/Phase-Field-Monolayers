#Run this either from the folder:
#-- of a simulation run's data, passing NO EXTRA ARGUMENTS
#-- from which the simulation was called, passing the data directory's name as the only argument
#runner automatically calls this appropriately

#Uses "Pillow": you may need to "python3 -m pip install --upgrade Pillow" (or just "python", instead of "python3")
#More info on https://pillow.readthedocs.io/en/stable/installation.html
#Note the WARNING: "Pillow and PIL cannot co-exist in the same environment"

import glob
from PIL import Image
import sys

gifName = 'out.gif'

gifFrameDuration = 20
stride = 10
initialFramesWithoutJumping = 10

def pgmToColoredJpg(source_path, dest_path):
    im = Image.open(source_path)
    im = im.convert("RGB")

    source = im.split()
    R, G, B = 0, 1, 2

    fullBlue = 0
    fullBlack = 18
    fullWhite = 237
    fullRed = 255

    maskFullBlue = source[R].point(lambda i: i <= fullBlue and 255)
    maskBlueBlack = source[R].point(lambda i: i > fullBlue and i < fullBlack and 255)
    maskBlackWhite = source[R].point(lambda i: i >= fullBlack and i < fullWhite and 255)
    maskWhiteRed = source[R].point(lambda i: i >= fullWhite and i < fullRed and 255)
    maskFullRed = source[R].point(lambda i: i >= fullRed and 255)

    newBlue = source[B].point(lambda i: 255)

    newBlueBlack = source[B].point(lambda i: 255*(fullBlack - i)/(fullBlack - fullBlue))
    
    newBlackWhite = source[G].point(lambda i: 255*(i - fullBlack)/(fullWhite - fullBlack))    

    newWhiteNonRed = source[R].point(lambda i: 255*(1 - ((i - fullWhite)/(fullRed - fullWhite))))

    newRed = source[R].point(lambda i: 255)

    source[R].paste(newBlue, None, maskFullBlue)
    source[G].paste(0, None, maskFullBlue)
    source[B].paste(0, None, maskFullBlue)

    source[R].paste(0, None, maskBlueBlack)
    source[G].paste(0, None, maskBlueBlack)
    source[B].paste(newBlueBlack, None, maskBlueBlack)

    source[R].paste(newBlackWhite, None, maskBlackWhite)
    source[G].paste(newBlackWhite, None, maskBlackWhite)
    source[B].paste(newBlackWhite, None, maskBlackWhite)

    source[R].paste(255, None, maskWhiteRed)
    source[G].paste(newWhiteNonRed, None, maskWhiteRed)
    source[B].paste(newWhiteNonRed, None, maskWhiteRed)

    source[R].paste(newRed, None, maskFullRed)
    source[G].paste(0, None, maskFullRed)
    source[B].paste(0, None, maskFullRed)

    imNew = Image.merge(im.mode, source)

    imNew.save(dest_path, "JPEG", optimize=True, quality=100)

if len(sys.argv) > 2:
    print('received bad arguments')
    exit
    
dataPath = ''
if len(sys.argv) == 2:
    dataPath = sys.argv[1] + '\\'

pgmPaths = glob.glob(dataPath + "*.pgm")
for path in pgmPaths:
    pgmToColoredJpg(path, path[:-4] + ".jpg")
    
jpgPaths = sorted(glob.glob(dataPath + "*.jpg"), key = len)
images = []
for path in jpgPaths:
    if frame < initialFramesWithoutJumping or (frame % stride) == 0:
        im = Image.open(path)
        images.append(im)
    frame += 1
images[0].save(dataPath + gifName, save_all=True, append_images=images[1:], 
               optimize=True, duration=gifFrameDuration, loop=0, comment=jpgPaths[0])

print('=> jpgs and gif available on the data folder (' + dataPath + ')\n')