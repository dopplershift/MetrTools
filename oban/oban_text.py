import ob
import numpy as N

def get_obs(filename):
  '''Reads in one of the retardedly formatted text files where information for
     a station spans more than one line.'''
  obInfo = []
  for line in open(filename, 'r'):
    obInfo.extend(line.split())
  obs = []
  for index in range(0,len(obInfo),6):
    data = ob.ob()
    data.stid = float(obInfo[index])
    data.lat = float(obInfo[index+1])
    data.lon = float(obInfo[index+2])
    data.height = float(obInfo[index+3])
    data.dir = float(obInfo[index+4])
    data.speed = float(obInfo[index+5])
    obs.append(data)
  return obs

def get_arrays(filename):
  '''Reads in one of the retardedly formatted text files where information for
     a station spans more than one line. Data is returned as a collection of
     arrays'''
  class collection(object): pass
  data = collection()
  data_table = N.fromfile(filename, sep=' ').reshape(-1,6)
  data.stids = data_table[:,0]
  data.lats = data_table[:,1]
  data.lons = data_table[:,2]
  data.heights = data_table[:,3]
  data.wind_dirs = data_table[:,4]
  data.wind_speeds = data_table[:,5]
  return data
