import sys
import os

fxn = sys.argv[1]
[module,method] = fxn.split(".")
module = __import__("%s" %module)

# Run the command
args = ["sys.argv["+str(i)+"]" for i in range(2,len(sys.argv))]
eval("module."+method+"("+",".join(args)+")")
