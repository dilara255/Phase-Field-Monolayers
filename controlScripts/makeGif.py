#Run this either from the folder:
#-- of a simulation run's data, passing NO EXTRA ARGUMENTS
#-- from which the simulation was called, passing the data directory's name as the only argument
#runner automatically calls this appropriately

#Uses "Pillow": you may need to "python3 -m pip install --upgrade Pillow"
#More info on https://pillow.readthedocs.io/en/stable/installation.html
#Note the WARNING: "Pillow and PIL cannot co-exist in the same environment"

import glob
from PIL import Image
import sys

gifName = 'out.gif'
gifFrameDuration = 40

def pgmToColoredJpg(source_path, dest_path):
    im = Image.open(source_path)
    im = im.convert("RGB")

    source = im.split()
    R, G, B = 0, 1, 2

    fullBlack = 0
    fullBlue = 66
    fullYellow = 128
    fullRed = 190
    fullWhite = 255

    maskFullBlack = source[R].point(lambda i: i <= fullBlack and 255)
    maskBlackBlue = source[R].point(lambda i: i > fullBlack and i < fullBlue and 255)
    maskBlueYellow = source[R].point(lambda i: i >= fullBlue and i < fullYellow and 255)
    maskYellowRed = source[R].point(lambda i: i >= fullYellow and i < fullRed and 255)
    maskRedWhite = source[R].point(lambda i: i >= fullRed and i < fullWhite and 255)
    maskFullWhite = source[R].point(lambda i: i >= fullWhite and 255)

    newBlack = source[R].point(lambda i: 0)

    newBlackBlue = source[B].point(lambda i: 255*(i - fullBlack)/(fullBlue - fullBlack))

    newBlueYellowBlue = source[B].point(lambda i: 255*(1-(i - fullBlue)/(fullYellow - fullBlue)))
    newBlueYellowYellow = source[B].point(lambda i: 255*(i - fullBlue)/(fullYellow - fullBlue))

    newYellowRedGreen = source[B].point(lambda i: 255*(fullRed - i)/(fullRed - fullYellow))

    newRedWhiteGreenBlue = source[R].point(lambda i: 255*(i - fullRed)/(fullWhite - fullRed))

    newWhite = source[R].point(lambda i: 255)

    source[R].paste(newBlack, None, maskFullBlack)
    source[G].paste(newBlack, None, maskFullBlack)
    source[B].paste(newBlack, None, maskFullBlack)

    source[R].paste(newBlack, None, maskBlackBlue)
    source[G].paste(newBlack, None, maskBlackBlue)
    source[B].paste(newBlackBlue, None, maskBlackBlue)

    source[R].paste(newBlueYellowYellow, None, maskBlueYellow)
    source[G].paste(newBlueYellowYellow, None, maskBlueYellow)
    source[B].paste(newBlueYellowBlue, None, maskBlueYellow)

    source[R].paste(newWhite, None, maskYellowRed)
    source[G].paste(newYellowRedGreen, None, maskYellowRed)
    source[B].paste(newBlack, None, maskYellowRed)

    source[R].paste(newWhite, None, maskRedWhite)
    source[G].paste(newRedWhiteGreenBlue, None, maskRedWhite)
    source[B].paste(newRedWhiteGreenBlue, None, maskRedWhite)

    source[R].paste(newWhite, None, maskFullWhite)
    source[G].paste(newWhite, None, maskFullWhite)
    source[B].paste(newWhite, None, maskFullWhite)

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
    im = Image.open(path)
    images.append(im)
images[0].save(dataPath + gifName, save_all=True, append_images=images[1:], 
               optimize=True, duration=gifFrameDuration, loop=0, comment=jpgPaths[0])

print('=> jpgs and gif available on the data folder (' + dataPath + ')\n')