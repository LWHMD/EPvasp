#!/usr/bin/python3
# -*- coding: UTF-8 -*-

########### script to extract data from PROCAR ############
# Input file : PROCAR, 278                                #
#              KPOINTS,                                   #
#              POSCAR, fermi.dat                          #
########### script to extract data from EIGENVAL ##########
# Input file : EIGENVAL,(soc,nosoc.megnetic)              #
#              KPOINTS,                                   #
#              POSCAR, fermi.dat                          #
#              KPOINTS.DFT(HSE),                          #
# run command: python3.4 VaspToolX.py                     #
# Author     : Leiwang  updata 2019/07/18                #
###########################################################
# The version copy from ubuntu

# note that the format of fermi.dat
#   ISMEAR =     0;   SIGMA  =   0.01  broadening in eV -4-tet -1-fermi 0-gaus
# E-fermi :   7.0717     XC(G=0): -11.2821     alpha+bet :-12.3742

# creat it by : grep fermi OUTCAR > fermi.dat

#import numpy as np
#$from sympy import *
import math
import os

pi = math.pi
sqrt = math.sqrt
#  used to read data from file
def read_data(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
    return content

#  used to write data to file
def write2txt(filename, data):
    f = open(filename, 'a')
    f.write(data + "\n")
    f.close()

#read fermi level
def fermienergy(fermi):
    efermi_tem = read_data(fermi)
    efermi0 = efermi_tem[1]
    efermi = efermi0.split()
    efermi_final = efermi[2]
    return efermi_final

# tranfer real space lattice a1 a2 a3 to reciprocal space lattice b1 b2 b3
def real_to_reciprocal(lines0):

    a1 = lines0[2]
    a2 = lines0[3]
    a3 = lines0[4]
    a1 = a1.split()
    a2 = a2.split()
    a3 = a3.split()
    #print(len(a1))
    #print(a1[0])
    A = [[],[],[]]
    A[0].append(a1[0])
    A[0].append(a1[1])
    A[0].append(a1[2])
    A[1].append(a2[0])
    A[1].append(a2[1])
    A[1].append(a2[2])
    A[2].append(a3[0])
    A[2].append(a3[1])
    A[2].append(a3[2])
    #print(A)
    volume = (float(a1[0])*float(a2[1])*float(a3[2])+float(a1[1])*float(a2[2])*float(a3[0])
              +float(a1[2])*float(a2[0])*float(a3[1])-float(a1[0])*float(a2[2])*float(a3[1])
              -float(a1[1])*float(a2[0])*float(a3[2])-float(a1[2])*float(a2[1])*float(a3[0]))
    b=[[],[],[]]
    c=[]
    for i in (0,1,2):
        if i==0:
            j = 1
            k = 2
        elif i==1:
            j = 2
            k = 0
        else:
            j = 0
            k = 1
        c.append(float(A[j][1])*float(A[k][2])-float(A[j][2])*float(A[k][1]))
        c.append(float(A[j][2])*float(A[k][0])-float(A[j][0])*float(A[k][2]))
        c.append(float(A[j][0])*float(A[k][1])-float(A[j][1])*float(A[k][0]))
        #print (c)
        for l in (0,1,2):
            bx = 2*pi*float(c[l])/volume
            b[i].append(bx)
        #print(b[i])
        del c[:]
    return b

# calculate the distance between two point in k space
def L_in_kspace(ary1,ary2,b):
    # get first point
    b1=[[],[],[]]
    b2=[[],[],[]]
    b3=[[],[],[]]
    b1[0].append(b[0][0]*ary1[0])
    b1[0].append(b[0][1]*ary1[0])
    b1[0].append(b[0][2]*ary1[0])
    b1[1].append(b[1][0]*ary1[1])
    b1[1].append(b[1][1]*ary1[1])
    b1[1].append(b[1][2]*ary1[1])
    b1[2].append(b[2][0]*ary1[2])
    b1[2].append(b[2][1]*ary1[2])
    b1[2].append(b[2][2]*ary1[2])
    #print (b1)
    # get second point
    b2[0].append(b[0][0]*ary2[0])
    b2[0].append(b[0][1]*ary2[0])
    b2[0].append(b[0][2]*ary2[0])
    b2[1].append(b[1][0]*ary2[1])
    b2[1].append(b[1][1]*ary2[1])
    b2[1].append(b[1][2]*ary2[1])
    b2[2].append(b[2][0]*ary2[2])
    b2[2].append(b[2][1]*ary2[2])
    b2[2].append(b[2][2]*ary2[2])
    #print (b2)
    # get the vector first point to second point
    b3[0].append(b2[0][0]-b1[0][0])
    b3[0].append(b2[0][1]-b1[0][1])
    b3[0].append(b2[0][2]-b1[0][2])
    b3[1].append(b2[1][0]-b1[1][0])
    b3[1].append(b2[1][1]-b1[1][1])
    b3[1].append(b2[1][2]-b1[1][2])
    b3[2].append(b2[2][0]-b1[2][0])
    b3[2].append(b2[2][1]-b1[2][1])
    b3[2].append(b2[2][2]-b1[2][2])
    #print(b3)

    # to get the mod of vector
    kb1 = 1/(2*pi)*sqrt((b3[0][0]+b3[1][0]+b3[2][0])**2+(b3[0][1]+b3[1][1]+b3[2][1])**2+(b3[0][2]+b3[1][2]+b3[2][2])**2)
    return kb1      #

# To calculate k mesh
def calcu_k_meth(lines0,lines1):
    result = []
    for line in lines1:                          #to read each line
        line = line.strip()                             #去掉每行头尾空白
        if not len(line) or line.startswith('#'):       #判断是否是空行或注释行
             continue                                    #是的话，跳过不处理
        result.append(line)                             #保存
    mesh = int(result[1])
    #print(mesh)
    i = 4     # initial line
    K_path = []
    while i < len(result):
        K_path.append(result[i])
        i += 1

    # get mesh
    Nk_path = len(K_path)
    #print(Nk_path)
    L_k_tem = 0.0
    L_k_mesh_list = []
    h_k = []    # high symmetry point
    for j in range(0,Nk_path,2):
        p1 = K_path[j]
        p1 = p1.split()
        p3=[]
        for char in p1:
            char = float(char)
            p3.append(char)
        p2 = K_path[j+1]
        p2 = p2.split()
        p4=[]
        for char in p2:
            char = float(char)
            p4.append(char)
        #print(p3,p4)
        reci = real_to_reciprocal(lines0)
        #print (reci)
        L_k = L_in_kspace(p3,p4,reci)     # calculate the distance between two point in k space
        for i in range(0,mesh,1):
            L_k_mesh = (L_k)*i/(mesh-1)
            L_k_mesh_list.append(L_k_mesh+L_k_tem)
        h_k.append(L_k_tem)
        L_k_tem = L_k_tem+L_k
    returnterm=[]
    returnterm.append(L_k_mesh_list)
    returnterm.append(mesh)
    return returnterm

# used to calculate high symmetry line
def high_symmetry_line(lines0,lines1):
    k_mesh_reci0=calcu_k_meth(lines0,lines1)
    k_mesh_reci1=k_mesh_reci0[0]
    k_mesh = k_mesh_reci0[1]
    kpoint_high_sym=[]
    i=0
    kpoint_high_sym.append(k_mesh_reci1[i])
    while i <len(k_mesh_reci1):
        i=i+k_mesh
        kpoint_high_sym.append(k_mesh_reci1[i-1])
    return kpoint_high_sym

# Deal with PROCAR file and get the orbital component
def project_orbit():

    while True:
        conform_file = str(input('To ensure POSCAR, PROCAR, KPOINTS, fermi.dat in current floder: Y/N'))
        if  'Y' != conform_file :
            print('please prepare POSCAR, PROCAR, KPOINTS ')
            continue
        else:
            break


    lines0 = read_data('POSCAR')     #read  POSCAR
    lines1 = read_data('KPOINTS')
    lines3 = read_data('PROCAR')

    # begin to extract the orbit component from PROCAR
    #print (lines0[5])
    #while True:
    #element = str(input('input the kind of element:'))
    #    if element not in lines0[5]:
    #        print('the element is not right, please input it again')
    #        continue
    #    else:
    #        break



    # extract data in two mode soc or nosoc
    #mode = int(input('spd input 1; s px py pz dxy dyz dz2 dxz dx2 input 2:'))   # LORBIT
    mode = 2
    LSO = int(input('nosoc input 1; soc input 2:'))     # LSORBIT
    mag = int(input('nonmagnetic 1; magnetic 2:'))   # ISPIN equal to 1 or 2
    efermi = fermienergy('fermi.dat')
    efermi = float(efermi)
    #nbands = int(input('the index of band : '))
    #filename = 'PROCAR.txt'  # both .txt file and this script at same content, in this way you don't need  to write the path
    #print (lines[0],'\t',lines[1],'\t',lines[2],'\t',lines[3])
    lines2 = lines3[1]
    lines2 = lines2.split()
    nk = int(lines2[3])    # read the number of kpoints from PROCAR
    nb = int(lines2[7])    # read the number of bands from PROCAR
    ni = int(lines2[11])    # read the number of ion from PROCAR
    #print(len(lines))
    print ('number of kpoints:',nk,'number of bands:',nb,'number of ion:',ni)

    L_k_mesh_list = calcu_k_meth(lines0,lines1)
    L_k_mesh_list=L_k_mesh_list[0]

    if LSO ==1:
        tb_betw=(ni+1)+4         # the number of line between two adjacent band in one k-block
    else:
        tb_betw=(ni+1)*4+4

    N_A = 0
    N_i = 0
    Num_A = []
    #if element in lines0[5]:
    Element=lines0[5].split()
    Num_A=lines0[6].split()
        #print (Element)
        #print (Num_A)
    i = 0
    while i< len(Element):
        N_A = N_A + int(Num_A[i])
        for m in range(0,mag,1):
            #print('m',m)
            for i_nb in range(0,nb,1):  #bands
                for i_nk in range(0,nk,1):   #kpoints
                #print(i)
                    nkblock = tb_betw*nb+3    # the number of line between two adjacent k-block, such k-points 1 and k-points 2
                    #print('nkblock:',nkblock)
                    k_tmp = lines3[3+i_nk*nkblock]        # the  fractional coordinate of k-points
                    k = k_tmp[19:52]
                    A = N_A-int(Num_A[N_i])+1
                    s = 0;p = 0;d = 0
                    px=0;py=0;pz=0;dxy=0;dyz=0;dxz=0;dx2=0;dz2=0
                    Energy = lines3[i_nk*nkblock+2+(tb_betw*(i_nb)+3)+m*(nk*nkblock+1)]
                    #print(Energy)
                    Energy = Energy.split()
                    energy = float(Energy[4])-efermi
                    #print(Energy)
                    if mode == 1:
                        for j in range(A,N_A+1,1):          #  To choose the line the atom that you choose located in
                            #print (j)
                            xx_tmp = lines3[i_nk*nkblock+2+(tb_betw*(i_nb)+3)+j+2+m*(nk*nkblock+1)]      # the line include the atom that you choose under nk,nb
                            #print (xx_tmp)
                            xx = xx_tmp.split()
                            s = s + float(xx[1])       #    s
                            #p = p + float(xx[2])+float(xx[3])+float(xx[4])   #     py     pz     px
                            p = p + float(xx[2])+float(xx[3])+float(xx[4])
                            d = d + float(xx[5])+ float(xx[6])+float(xx[7])+float(xx[8])+float(xx[9])   #    dxy    dyz    dz2    dxz    dx2
                            #d = d + float(xx[5])+ float(xx[6])+float(xx[7])+float(xx[8])+float(xx[9])
                        #write2txt('band-spd-'+element+'.txt',str(i_nb+1)+'\t'+str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(s)+'\t'+str(p)+'\t'+str(d))
                        write2txt('band-s-'+Element[i]+'.dat',str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(s))
                        write2txt('band-pxpy-'+Element[i]+'.dat',str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(p))
                        write2txt('band-d-'+Element[i]+'.dat',str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(d))
                    else:
                        for j in range(A,N_A+1,1):
                            #print (j)
                            xx_tmp = lines3[i_nk*nkblock+5+tb_betw*(i_nb)+j+2+m*(nk*nkblock+1)]
                            #print (xx_tmp)
                            xx = xx_tmp.split()
                            s = s + float(xx[1])
                            px = px + float(xx[2])
                            py = py + float(xx[3])
                            pz = pz + float(xx[4])
                            dxy = dxy + float(xx[5])
                            dyz = dyz + float(xx[6])
                            dz2 = dz2 +float(xx[7])
                            dxz = dxz +float(xx[8])
                            dx2 = dx2 +float(xx[9])
                        write2txt('band-spxdx-'+Element[i]+'.dat',str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(s)+'\t'+str(px)+'\t'+str(py)+'\t'+str(pz)+'\t'+str(dxy)+'\t'+str(dyz)+'\t'+str(dz2)+'\t'+str(dxz)+'\t'+str(dx2))
                        #write2txt('band-spxdx-'+element+'.dat',str(i_nb+1)+'\t'+str(L_k_mesh_list[i_nk])+'\t'+str(energy)+'\t'+str(s)+'\t'+str(px)+'\t'+str(py)+'\t'+str(pz)+'\t'+str(dxy)+'\t'+str(dyz)+'\t'+str(dz2)+'\t'+str(dxz)+'\t'+str(dx2))
                if mode == 1:
                    #write2txt('band-spd-'+element+'.txt',str( )+'\t')    # space
                    write2txt('band-s-'+Element[i]+'.dat',str( )+'\t')
                    write2txt('band-pxpy-'+Element[i]+'.dat',str( )+'\t')
                    write2txt('band-d-'+Element[i]+'.dat',str( )+'\t')
                else:
                    write2txt('band-spxdx-'+Element[i]+'.dat',str( )+'\t')

        i += 1
    hsl=high_symmetry_line(lines0,lines1)
    for i in range(len(hsl)):
        write2txt('high-symmetry-line.dat',str(hsl[i])+'\t'+str(-30))
        write2txt('high-symmetry-line.dat',str(hsl[i])+'\t'+str(30))
        write2txt('high-symmetry-line.dat',' ')
    write2txt('high-symmetry-line.dat',str(0)+'\t'+str(0))
    write2txt('high-symmetry-line.dat',str(hsl[len(hsl)-1])+'\t'+str(hsl[0]))

    project_orbit2()


def project_orbit2():
    print('This part is used to operate the orbit data in PROCAR')
    print('To choose the element')
    structure = read_data('POSCAR')
    print (structure[5])
    element0 = str(input('input the kind of element:'))
    element = element0.split()
    i = 0
    Norbit = []
    while i < len(element):
        print ('1. s 2. py 3. pz 4.px 5. dxy 6. dyz 7.dz2 8. dxz 9. x2-y2')
        Norbit0 = str(input('input the orbit of element'+str(element[i])+ """in format '1 2 3 4'"""))
        Norbit = Norbit0.split()
        write2txt('projected_band.dat','1. s 2. py 3. pz 4.px 5. dxy 6. dyz 7.dz2 8. dxz 9. x2-y2')
        write2txt('projected_band.dat','element : '+str(element[i])+'\t'+Norbit0)
        i += 1
    N_el = 0
    #print ('Nor',len(Norbit))
    while N_el < len(element):
        orbit_file0 = read_data('band-spxdx-'+element[N_el]+'.dat')
        for orbit_file in  orbit_file0:
            orbit = orbit_file.split()
            if len(orbit)==0:
                write2txt('projected_band.dat',str('')+'\t')
                continue
            path = orbit[0]
            energy = orbit[1]
            component = 0
            i=0

            while i < len(Norbit):
                N = int(Norbit[i])+1
                component = component + float(orbit[N])
                i += 1
            write2txt('projected_band.dat',str(path)+'\t'+str(energy)+'\t'+str(component)+'\t')
        N_el += 1





# used  to read EIGENVAL file
def read_eigenval(lines3,nk,nb,mag):
    eigenval_file = lines3
    #print (eigenval_file[7])
    i= 7
    list_eigenval_total=[[0 for i in range(nk)] for j in range(nb)]
    list_eigenval_up=[[0 for i in range(nk)] for j in range(nb)]
    list_eigenval_down=[[0 for i in range(nk)] for j in range(nb)]
    k=0
    if mag ==1:
        while i < len(eigenval_file):
            i = i + 1  # add one line
            for j in range(0,nb,1):
                value = eigenval_file[i]
                temp = value.split()
                i +=1
                if k==nk:
                    k=0
                list_eigenval_total[j][k]=temp[1]
            k+=1    # index of k points
            i+=1    # add one line
        return list_eigenval_total
    else:
        while i < len(eigenval_file):
            i = i + 1  # add one line
            for j in range(0,nb,1):
                value = eigenval_file[i]
                temp = value.split()
                i +=1
                if k==nk:
                    k=0
                list_eigenval_up[j][k]=temp[1]
            k+=1    # index of k points
            i+=1    # add one line
        i=7
        k=0
        while i < len(eigenval_file):
            i = i + 1  # add one line
            for j in range(0,nb,1):
                value = eigenval_file[i]
                temp = value.split()
                i +=1
                if k==nk:
                    k=0
                list_eigenval_down[j][k]=temp[2]
            k+=1    # index of k points
            i+=1    # add one line
        list = [list_eigenval_up,list_eigenval_down]
        return list

# used to calculate normal band structure
def band_cal():
    lines0 = read_data('POSCAR')     #read  POSCAR
    lines1 = read_data('KPOINTS')
    lines3 = read_data('EIGENVAL')
    mag = int(input('nonmagnetic 1; magnetic (nosoc) 2:'))   # ISPIN equal to 1 or 2
    efermi = fermienergy('fermi.dat')
    efermi = float(efermi)
    lines2 = lines3[5]
    lines2 = lines2.split()
    num_k = int(lines2[1])    # read the number of kpoints from PROCAR
    num_b = int(lines2[2])    # read the number of bands from PROCAR
    print ('number of kpoints:',num_k,'number of bands:',num_b)
        # extract data in two mode magnetic or no
    if 1 == mag:
        L_k_mesh_list = calcu_k_meth(lines0,lines1)
        L_k_mesh_list=L_k_mesh_list[0]
        list_eigen_val = read_eigenval(lines3,num_k,num_b,mag)
        for ib in range(num_b):
            for ik in range(num_k):
                write2txt('bandstructure.dat',str(L_k_mesh_list[ik])+'\t'+str(float(list_eigen_val[ib][ik])-efermi))
            write2txt('bandstructure.dat',' ')
        hsl=high_symmetry_line(lines0,lines1)
        for i in range(len(hsl)):   # print High symmetry line
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(-30))
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(30))
            write2txt('bandstructure.dat',' ')
        write2txt('bandstructure.dat',str(0)+'\t'+str(0))
        write2txt('bandstructure.dat',str(hsl[len(hsl)-1])+'\t'+str(hsl[0]))
    elif 2==mag :
        L_k_mesh_list = calcu_k_meth(lines0,lines1)
        L_k_mesh_list=L_k_mesh_list[0]
        list_eigen_val = read_eigenval(lines3,num_k,num_b,mag)
        list_eigen_val_up = list_eigen_val[0]
        list_eigen_val_down = list_eigen_val[1]

        for ib in range(num_b):
            for ik in range(num_k):
                write2txt('bandstructure.dat',str(L_k_mesh_list[ik])+'\t'+str(float(list_eigen_val_up[ib][ik])-efermi))
            write2txt('bandstructure.dat',' ')
        for ib in range(num_b):
            for ik in range(num_k):
                write2txt('bandstructure.dat',str(L_k_mesh_list[ik])+'\t'+str(float(list_eigen_val_down[ib][ik])-efermi))
            write2txt('bandstructure.dat',' ')
        hsl=high_symmetry_line(lines0,lines1)
        for i in range(len(hsl)):   # print High symmetry line
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(-30))
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(30))
            write2txt('bandstructure.dat',' ')
        write2txt('bandstructure.dat',str(0)+'\t'+str(0))
        write2txt('bandstructure.dat',str(hsl[len(hsl)-1])+'\t'+str(hsl[0]))
    else:
        print('incorrect mogmam')

# used to read the kpoints of HSE calculation
def read_hse_KPOINTS(lines1,lines_1):
    result = []
    for line in lines_1:                          #to read each line
        line = line.strip()                             #去掉每行头尾空白
        if not len(line) or line.startswith('#'):       #判断是否是空行或注释行
             continue                                    #是的话，跳过不处理
        result.append(line)                             #保存
    mesh = int(result[1])
    #print(mesh)
    i = 4     # initial line
    K_path = []
    while i < len(result):
        K_path.append(result[i])
        i += 1

    # get mesh
    Nk_path = len(K_path)
    list = []     # used to store the point on the high symmetry line
    for j in range(0,Nk_path,2):
        p1 = K_path[j]
        p1 = p1.split()
        p3=[]
        for char in p1:
            char = float(char)
            p3.append(char)
        p2 = K_path[j+1]
        p2 = p2.split()
        p4=[]
        for char in p2:
            char = float(char)
            p4.append(char)
        #print(p3,p4)
        #print (reci)
        # output k points
        for i in range(mesh):
            list_k = []
            px = p3[0]-(p3[0]-p4[0])*(i)/(mesh-1)
            py = p3[1]-(p3[1]-p4[2])*(i)/(mesh-1)
            pz = p3[2]-(p3[2]-p4[2])*(i)/(mesh-1)
            list_k.append(px)
            list_k.append(py)
            list_k.append(pz)
            #print (list_k)
            list.append(list_k)
    #print (list)

    # compare with HSE mesh
    k_mesh_hse = []
    for i in range(3,len(lines1),1):
        kp0=lines1[i]
        kp1 = kp0.split()
        if kp1[3] == '0':
            #print (kp1)
            k_mesh_hse.append(kp1)
    #print (len(list),len(k_mesh_hse),k_mesh_hse)
    #print(list)
    com_num = []     # used to collect compare number
    list_temp = list
    k=0
    j=0
    for i in range(len(k_mesh_hse)):
        #print ('i',i)
        while j<len(list):
            com_hse_kx = k_mesh_hse[i][0]
            com_hse_ky = k_mesh_hse[i][1]
            com_dft_kx = list_temp[j][0]
            com_dft_ky = list_temp[j][1]
            #print(k_mesh_hse[i],list_temp[j])
            if abs(float(com_hse_kx)-com_dft_kx) <= 0.00001 and abs(float(com_hse_ky)-com_dft_ky <= 0.00001):
                com_num.append(j)
                k=j
                #print ('k',k)
                break
            j += 1
        j=k+1
    #print(len(com_num),com_num)
    return com_num

# used to read EIGENVAL of HSE calculation
def read_hse_eigenval(lines3,nk,nb):
    eigenval_file = lines3
    i= len(eigenval_file)-nk*(nb+2)+1
    list_eigenval_total=[[0 for i in range(nk)] for j in range(nb)]
    k=0
    #print (i,eigenval_file[i])
    while i < len(eigenval_file):
        i = i + 1  # add one line
        for j in range(0,nb,1):
            value = eigenval_file[i]
            temp = value.split()
            i +=1
            if k==nk:
                k=0
            list_eigenval_total[j][k]=temp[1]
        k+=1    # index of k points
        i+=1    # add one line
    return list_eigenval_total

# used to calculate hse band
def band_hse_cal():
    lines0 = read_data('POSCAR')     #read  POSCAR
    lines1 = read_data('KPOINTS')
    lines3 = read_data('EIGENVAL')
    print('be sure set correct number of kpoints in EIGENVAL file')
    while True:
        k_DFT = input('input KPOINTS.DFT in current floder (Y/N):')
        if k_DFT != 'Y':
            print('please input normal mode KPOINTS.DFT file')
        else:
            print ('OK')
            break

    lines_1 = read_data('KPOINTS.DFT')
    print('1. normal mode input')
    print('2. abnormal mode input')
    mode = float(input())
    efermi = fermienergy('fermi.dat')
    efermi = float(efermi)

    lines2 = lines3[5]    # to get the number of band from EGIENVAL file
    lines2 = lines2.split()
    num_b = int(lines2[2])

    k_mesh_hse_num = []
    for i in range(3,len(lines1),1):
        kp0=lines1[i]
        kp1 = kp0.split()
        if kp1[3] == '0':
            #print (kp1)
            k_mesh_hse_num.append(kp1)
    num_k = len(k_mesh_hse_num) # read the number of kpoints from EIGENVAL
    print(num_k)
    print('number of kpoints:',num_k)
    print('number of bands:',num_b)
    # extract data in two mode magnetic or no

    if mode == 1:
        L_k_mesh_list = calcu_k_meth(lines0,lines_1)
        L_k_mesh_list=L_k_mesh_list[0]
        list_eigen_val = read_hse_eigenval(lines3,num_k,num_b)

        for ib in range(num_b):
            for ik in range(num_k):
                write2txt('bandstructure.dat',str(L_k_mesh_list[ik])+'\t'+str(float(list_eigen_val[ib][ik])-efermi))
            write2txt('bandstructure.dat',' ')
        hsl=high_symmetry_line(lines0,lines_1)
        for i in range(len(hsl)):   # print High symmetry line
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(-30))
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(30))
            write2txt('bandstructure.dat',' ')
        write2txt('bandstructure.dat',str(0)+'\t'+str(0))
        write2txt('bandstructure.dat',str(hsl[len(hsl)-1])+'\t'+str(hsl[0]))
    else:
        L_k_mesh_list = calcu_k_meth(lines0,lines_1)
        L_k_mesh_list = L_k_mesh_list[0]
        L_k_mesh_list_com = []
        num = read_hse_KPOINTS(lines1,lines_1)
        #print(len(num))
        for i in range(len(num)):
            #print(i,int(num[i]))
            L_k_mesh_list_com.append(L_k_mesh_list[int(num[i])])
        list_eigen_val = read_hse_eigenval(lines3,num_k,num_b)
        for ib in range(num_b):
            for ik in range(len(num)):
                write2txt('bandstructure.dat',str(L_k_mesh_list_com[ik])+'\t'+str(float(list_eigen_val[ib][ik])-efermi))
            write2txt('bandstructure.dat',' ')
        hsl=high_symmetry_line(lines0,lines_1)
        for i in range(len(hsl)):   # print High symmetry line
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(-30))
            write2txt('bandstructure.dat',str(hsl[i])+'\t'+str(30))
            write2txt('bandstructure.dat',' ')
        write2txt('bandstructure.dat',str(0)+'\t'+str(0))
        write2txt('bandstructure.dat',str(hsl[len(hsl)-1])+'\t'+str(hsl[0]))

# used to calculate band structure
def bandstructure():
    conform_file = str(input('To ensure POSCAR, EIGENVAL, KPOINTS, fermi.dat in current floder: Y/N'))
    if  'Y' == conform_file :
        print('please prepare POSCAR, EIGENVAL, KPOINTS ')
        print('To choose the program that you want to use: ')
        print('1. normal band')
        print('2. HSE band')
        choose_mode = str(input())
        if '1' ==choose_mode:
            band_cal()
        else:
            band_hse_cal()

def band_kpoint_PROCAR():
    LSO = int(input('nosoc 1 or soc 2'))
    ONE_kpoint = int(input('input one k-point'))
    SOME_bands0 = str(input('input bands'))
    SOME_bands=SOME_bands0.split()
    procar = read_data('PROCAR')
    procar_line2 = procar[1]
    kpoints_bands_ions = procar_line2.split()
    kpoints = int(kpoints_bands_ions[3])
    bands = int(kpoints_bands_ions[7])
    ions = int(kpoints_bands_ions[11])
    print ('number of kpoints:',kpoints,'number of bands:',bands)
    i=0
    j=0
    # To find the
    procar_line=''
    for procar_line in procar:
        procar_line_detail = procar_line.split()
        if 'k-point ' in procar_line and ONE_kpoint == int(procar_line_detail[1]) :
                #print (procar_line_detail[1])
                j=i
        i+=1
    kpoint_detail=[]
    block = 2+bands*(4+(ions+1)*(LSO**2))-1
    #print (j)
    for i in range(j-1,j+block-1,1):
        kpoint_detail.append(procar[i])
    write2txt('procar_bands_kpoint.dat','k-points :'+'\t'+str(ONE_kpoint))
    write2txt('procar_bands_kpoint.dat','bands :'+'\t'+str(SOME_bands0))
    ORBIT =procar[j+4]
    ORBIT = ORBIT[:-1]
    write2txt('procar_bands_kpoint.dat',ORBIT)
    i=0
    k=0
    #print (kpoint_detail)
    while i < len(SOME_bands):
        j=0
        for component_line in kpoint_detail:
            component=component_line.split()
            if 'band ' in component_line and str(SOME_bands[i]) == component[1]:
                k=j
            j+=1
            #print (j)
        i+=1
        print(k)
        bandsx=kpoint_detail[k+ions+3]
        bandsx=bandsx[:-1]
        write2txt('procar_bands_kpoint.dat',bandsx)


    # used to choose the mode you want to calculate
while True:
    print('To choose the program that you want to use:')
    print('1. project orbit')
    print('2. band structure')
    print('3. the component of some bands at one k-point')
    print('4. quit')
    project = str(input())
    if  '1' == project :
        print('you are performing a project-orbit program now.')
        project_orbit()
        continue
    elif  project == '2':
        bandstructure()
    elif  project == '3':
        band_kpoint_PROCAR()
    else:
        break
