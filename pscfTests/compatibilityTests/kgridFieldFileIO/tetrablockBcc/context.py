
import sys
import os

# add parent directory to path
#print("At start:\n\t{}\n\n".format(sys.path))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#print("After 1st insert:\n\t{}\n\n".format(sys.path))

# If parent directory does not have project, recurse up file structure
maxRecurse = 5 # move at most 5 directory levels up
try:
    import psoinverse
except(ImportError):
    notFound = True
    recurseLevel = 0
    while notFound and recurseLevel < maxRecurse:
        recurseLevel += 1
        notFound = False
        sys.path.insert(0, os.path.abspath(os.path.join(sys.path.pop(0), '..')))
        #print("At recurse level {}:\n\t{}\n\n".format(recurseLevel,sys.path))
        try:
            import psoinverse
        except(ImportError):
            notFound = True
            if recurseLevel >= maxRecurse:
                #sys.path.pop(0)
                #print("failed to find psoinverse",sys.path)
                raise(ImportError("Unable to find module psoinverse."))
        

# clean up python search path
#print("cleaning up path:\n\t{}\n\n".format(sys.path))
sys.path.pop(0)
#print("cleaned path:\n\t{}".format(sys.path))
