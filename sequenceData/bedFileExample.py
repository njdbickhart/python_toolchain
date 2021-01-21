import os
import sys

description = "A program to calculate the length and to then sort bed intervals."
usage = "{}\npython3 {} <Input bed file> <output sorted bed>".format(description, sys.argv[0])

def main(args):
    if len(args) != 2:
        # Sanity check to make sure we have a file
        print(usage)
        sys.exit()

    # We're going to define a list to contain all of the objects we want to store
    objList = []

    # Now to open our input and to create bedClass objects for each line of Input
    with open(args[0], 'r') as input:
        for l in input:
            s = l.rstrip().split()
            chr = s[0]
            start = s[1]
            end = s[2]
            # If there are additional columns in the bed file, keep them. Else, pass an empty list
            rest = s[3:]
            # This is our object, which is based on the "blueprint" of the class definition below
            # We have to pass at least three variables, but have room for more
            bedObject = bedClass(chr, start, end, rest)

            # Now we store the bedObject in our list of objects
            objList.append(bedObject)

    # Now, we're going to loop through the list and run all of our methods on each object
    for object in objList:
        # The method will work on the attributes in each object
        object.calcLen()
        # We are running this method without arguments becasue the default are OK
        object.createUCSC()

    # As promised, let's sort the list of bed objects
    # We ideally want to sort by chromosome and then by start coordinate
    # Both criteria need to be sorted at the same time, so we need to use a function that provides those keys
    # NOTE: we need to pass the actual function, so note that we don't include the parentheses at the end!
    # The bedSort funciton is defined below -- check for the "def bedSort"
    objList = sorted(objList, key=bedSort)

    # Our methods are finished, so let's print this out
    with open(args[1], 'w') as output:
        for object in objList:
            # Because we defined a method to create our output, it is super convenient to just print the output of that method!
            output.write(object.generateOutput())

    print("We're done here!")

def bedSort(bedObject):
    return bedObject.chr, bedObject.start

class bedClass:

    def __init__(self, chr, start, end, *args):
        # self should always be the first variable defined for all methods
        # Any other variables after "self" are initial variables that need to be present
        # to initialize the object. In this case, each bed object should at least have
        # a chr, start and end. The "args" are other elements of the bed entry that will be saved

        # Initial attributes
        self.chr = chr
        self.start = int(start)
        self.end = int(end)

        self.args = args if len(args) > 0 else []

        # this is optional, but I am predefining the attributes that we will fill later
        self.len = 0
        self.ucsc = ''

    def calcLen(self):
        # This is an example of a method that runs without any other Input
        # calculate the length of the bed interval
        self.len = self.end - self.start

    def createUCSC(self, fsep = ':', ssep = '-'):
        # And this is a method that can take alternative values for arguments
        # fsep is the first separator, and ssep is the second separator
        # This allows for customization of how to delimit the chromosome coordinate string
        # Again, "self" is required at first because that stores the attributes of the object!
        self.ucsc = f'{self.chr}{fsep}{self.start}{ssep}{self.end}'

    def generateOutput(self):
        # This is a convenience method to format all of the information we've collected
        # It creates a formatted string so that we can easily print out the class contents
        tempList = [self.chr, str(self.start), str(self.end)]
        # add any additional bed columns if present
        if not self.args:
            tempList.extend(self.args)
        tempList.append(str(self.len))
        tempList.append(str(self.ucsc))
        # Now to return our formatted string
        return "\t".join(tempList) + "\n"

if __name__ == "__main__":
    main(sys.argv[1:])
