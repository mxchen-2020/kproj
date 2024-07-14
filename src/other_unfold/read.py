import yaml

phonon_file=open('band.yaml','r',encoding='utf-8')
file_data = phonon_file.read()
phonon_file.close()
phonon_data=yaml.safe_load(file_data)
im=0
real=0
i=0

nqpoint=phonon_data['nqpoint']
natom=phonon_data['natom']

out_file=open("eigenvector.dat",'w')
while i<nqpoint:
    j=0
    while j<3*natom:
        k=0
        while k<natom:
            m=0
            while m<3:
                real= str(phonon_data['phonon'][i]['band'][j]['eigenvector'][k][m][0])
                im  = str(phonon_data['phonon'][i]['band'][j]['eigenvector'][k][m][1])
                print('(%s,%s)' % (real,im),file=out_file)
                m=m+1
            k=k+1
        j=j+1
    i=i+1
out_file.close()

kpoint_file=open("kpoint.dat",'w')

n=0

b11=phonon_data['reciprocal_lattice'][0][0]
b12=phonon_data['reciprocal_lattice'][0][1]
b13=phonon_data['reciprocal_lattice'][0][2]

b21=phonon_data['reciprocal_lattice'][1][0]
b22=phonon_data['reciprocal_lattice'][1][1]
b23=phonon_data['reciprocal_lattice'][1][2]

b31=phonon_data['reciprocal_lattice'][2][0]
b32=phonon_data['reciprocal_lattice'][2][1]
b33=phonon_data['reciprocal_lattice'][2][2]

while n<nqpoint:
    kpointx_D=0
    kpointy_D=0
    kpointz_D=0
    kpointx_D=(phonon_data['phonon'][n]['q-position'][0])
    kpointy_D=(phonon_data['phonon'][n]['q-position'][1])
    kpointz_D=(phonon_data['phonon'][n]['q-position'][2])
    kpointx_C= kpointx_D*b11+kpointy_D*b12+kpointz_D*b13
    kpointy_C= kpointx_D*b21+kpointy_D*b22+kpointz_D*b23
    kpointz_C= kpointx_D*b31+kpointy_D*b32+kpointz_D*b33
    print('%s %s %s' % (kpointx_C,kpointy_C,kpointz_C),file=kpoint_file )
    print('%s %s %s' % (kpointx_D,kpointy_D,kpointz_D),file=kpoint_file )
    n=n+1
kpoint_file.close()


distance_file=open("distance.dat","w")
n=0
while n<nqpoint:
    distance=phonon_data['phonon'][n]['distance']
    print(distance,file=distance_file)
    n=n+1
distance_file.close()

eigevalue_file=open('eigenvalue.dat','w')

i=0
while i<nqpoint:
    j=0
    while j<3*natom:
        E = phonon_data['phonon'][i]['band'][j]['frequency']
        print(E,file=eigevalue_file)
        j=j+1
    i=i+1




