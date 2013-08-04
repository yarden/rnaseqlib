import argh
from argh import arg

@arg("myarg", help="Argument")
def foo(myarg):
    print "Got passed %s" %(myarg)
