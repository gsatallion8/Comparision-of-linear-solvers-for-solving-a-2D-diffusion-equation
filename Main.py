import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from subprocess import call
from mpl_toolkits import mplot3d
import time

N=100
L=1.0
S=0.0
alpha=1.0
T1=400.0
T2=100.0
tol=0.1;

call(["gcc", "GS_method.c","-o","GS"])
call(["./GS",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol)])
call(["gcc", "lineTDMA_method.c","-o","lineTDMA"])
call(["./lineTDMA",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol)])
call(["gcc", "GD_method.c","-o","GD"])
call(["./GD",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol)])
call(["gcc", "CG_method.c","-o","CG"])
call(["./CG",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol)])

T=np.loadtxt("T_GS.txt")

x=np.linspace(0,L,N)
y=np.linspace(0,L,N)
fig=plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(x,y,T,100,cmap='binary')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('T')
fig.savefig("T_GS.png")

T=np.loadtxt("T_lineTDMA.txt")

x=np.linspace(0,L,N)
y=np.linspace(0,L,N)
fig=plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(x,y,T,100,cmap='binary')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('T')
fig.savefig("T_lineTDMA.png")

T=np.loadtxt("T_GD.txt")

x=np.linspace(0,L,N)
y=np.linspace(0,L,N)
fig=plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(x,y,T,100,cmap='binary')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('T')
fig.savefig("T_GD.png")

T=np.loadtxt("T_CG.txt")

x=np.linspace(0,L,N)
y=np.linspace(0,L,N)
fig=plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(x,y,T,100,cmap='binary')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('T')
fig.savefig("T_CG.png")

A=(T2-T1)/10.0
call(["gcc", "GS_method_IG.c","-lm","-o","GS_IG"])
call(["gcc", "lineTDMA_method_IG.c","-lm","-o","lineTDMA_IG"])
call(["gcc", "GD_method_IG.c","-lm","-o","GD_IG"])
call(["gcc", "CG_method_IG.c","-lm","-o","CG_IG"])
time_GS=[]
time_lineTDMA=[]
time_GD=[]
time_CG=[]
for k1  in range(1,5):
	time_GS.append([])
	time_lineTDMA.append([])
	time_GD.append([])
	time_CG.append([])
	for k2 in range(0,5):
		time_counter=time.time()
		call(["./GS_IG",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol),str(A),str(k1),str(k2)])
		time_GS[-1].append(time.time()-time_counter)
		time_counter=time.time()
		call(["./lineTDMA_IG",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol),str(A),str(k1),str(k2)])
		time_lineTDMA[-1].append(time.time()-time_counter)
		time_counter=time.time()
		call(["./GD_IG",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol),str(A),str(k1),str(k2)])
		time_GD[-1].append(time.time()-time_counter)
		time_counter=time.time()
		call(["./CG_IG",str(N),str(L),str(S),str(alpha),str(T1),str(T2),str(tol),str(A),str(k1),str(k2)])
		time_CG[-1].append(time.time()-time_counter)

		fig=plt.figure()
		err=np.loadtxt("GS_IG.txt")
		err.reshape((-1,2))
		plt.plot(err[:,0],err[:,1],color='r')
		err=np.loadtxt("lineTDMA_IG.txt")
		err.reshape((-1,2))
		plt.plot(err[:,0],err[:,1],color='b')
		err=np.loadtxt("GD_IG.txt")
		err.reshape((-1,2))
		if len(err)!=2:
			plt.plot(err[:,0],err[:,1],color='g')
		err=np.loadtxt("CG_IG.txt")
		err.reshape((-1,2))
		if len(err)!=2:
			plt.plot(err[:,0],err[:,1],color='y')
		string="Error_%d_%d.png" %(k1,k2)
		fig.savefig(string)
k2=np.array([0,1,2,3,4])
dictionary={0:"GS_time.png",1:"lineTDMA_time.png",2:"GD_time.png",3:"CG_time.png"}
k_label=["$k_1$=1","$k_1$=2","$k_1$=3","$k_1$=4"]
matrix=[time_GS,time_lineTDMA,time_GD,time_CG]
matrix=np.array(matrix)
for i in range(0,matrix.shape[0]):
	fig=plt.figure()
	for j in range(0,matrix.shape[1]):
		plt.plot(k2,matrix[i][j],label="%d" %j)
	plt.legend(k_label)
	fig.savefig(dictionary[i])
# np.savetxt("Times.txt",[time_GS,time_lineTDMA,time_GD,time_CG])
