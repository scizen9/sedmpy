import socket
import time
import os

from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_rawpath = parser.get('paths', 'rawpath')


class Reduce:
    """Class script to handle different commands run remotely from pharos"""
    def __init__(self):
        self.ip = "pharos.caltech.edu"
        self.port = 5006
        self.address = (self.ip,self.port)

    def sock_connect(self):
        s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        s.connect((self.address))

        return s

    def get_offset(self,abpair=False,file_name=""):
        """Get the a or ab offset from a file
        Format OFFSET,A or AB,Filename"""
        if abpair:
            abpair="AB"
        else:
            abpair="A"
        file_name = file_name.replace('s:/',_rawpath)
        cmd = "OFFSET,%s,%s\n" % (abpair,file_name)

        s = self.sock_connect()
        s.send(cmd)

        time.sleep(2)
        data = s.recv(1024)
        print "Got %s" % data

        s.close()
        return data

    def get_best_focus(self,files=""):
        """Send a list of focus images to get the best secondary focus
        position"""

        new_files = []
        for i in files:
            new_files.append(i.replace('s:/',_rawpath))

        cmd = "FOCUS,%s" % ",".join(new_files)

        s = self.sock_connect()
        s.send(cmd)

        time.sleep(1)
        data = s.recv(1024)
        print "Got %s" % data

        s.close()
        return data

    def get_sao_star(self,ra="",dec=""):
        "Get nearest SAO star"
        cmd ="SAO\n"

        s = self.sock_connect()
        s.send(cmd)

        time.sleep(1)
        data = s.recv(1024)
        print "Got %s" % data

        s.close()

        return data.split(',')

    def get_stats(self,file_name):
        file_name = file_name.replace('s:/',_rawpath)
        cmd = "STATS , %s\n" % file_name
        s = self.sock_connect()
        s.send(cmd)

        time.sleep(1)
        data = s.recv(1024)
        print "Got %s" % data

        s.close()


if __name__=='__main__':
	files = ['s:/20160305/rc20160305_03_15_40.fits',
	's:/20160305/rc20160305_03_16_21.fits',
	's:/20160305/rc20160305_03_17_02.fits',
	's:/20160305/rc20160305_03_17_43.fits',
	's:/20160305/rc20160305_03_18_24.fits',
	's:/20160305/rc20160305_03_19_05.fits',
	's:/20160305/rc20160305_03_19_46.fits']

	a = Reduce()
	#a.get_sao_star()
	a.get_offset(abpair=False,file_name='s:/20160304/rc20160304_07_42_07.fits')
	#a.get_best_focus(files)
	a.get_stats(file_name='s:/20160304/rc20160304_08_23_43.fits')


"""
s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
s.connect((ip,port))

s.send()

#files = glob.glob('s:/20160302/rc20160302_06_25_02.fits')
s.send('SOLVE %s\n' % '/20160302/rc20160302_06_25_02.fits' )
time.sleep(.3)
data = s.recv(2048)
print data
"""
