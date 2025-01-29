#Ultrafast Shape Recognition (USR) Algorithm
#------------------------------------------------------------------------------------------
import numpy as np
from ase.data import atomic_masses
from aegon.libutils import four_points
from aegon.libstin  import get_a_str, get_a_int, get_a_float
#------------------------------------------------------------------------------------------
log_file=get_a_str('output_file','glomos_out.txt')
verbose=get_a_int('verbose',0)
tole=get_a_float('tol_energy',0.1)
disc_file='discriminate_sym.xyz'
#------------------------------------------------------------------------------------------
def lastwo_central_moment(listxxx):
    xavg=np.mean(listxxx)
    res = listxxx - xavg
    ssq2 = np.mean(res**2) 
    ssq3 = np.mean(res**3) 
    return [ssq2,ssq3]
#------------------------------------------------------------------------------------------
def descriptors(moleculein):
####extra-descriptor USR algorithm: 
####standard USR algorithm +
####mass weighted USR algorithm
    listm=moleculein.get_masses()
    mavg = np.mean(listm)
    ctd, cst, fct, ftf=four_points(moleculein)
    lctd,  lcst,  lfct,  lftf = [], [], [], []
    lctdm, lcstm, lfctm, lftfm = [], [], [], []
    for iatom in moleculein:
        xyz=np.array(iatom.position)
        ictd1=np.linalg.norm(xyz-ctd)
        icst1=np.linalg.norm(xyz-cst)
        ifct1=np.linalg.norm(xyz-fct)
        iftf1=np.linalg.norm(xyz-ftf)
        lctd.append(ictd1)
        lcst.append(icst1)
        lfct.append(ifct1)
        lftf.append(iftf1)
        mi=atomic_masses[iatom.number]
        ictdm1=ictd1*mi/mavg
        icstm1=icst1*mi/mavg
        ifctm1=ifct1*mi/mavg
        iftfm1=iftf1*mi/mavg
        lctdm.append(ictdm1)
        lcstm.append(icstm1)
        lfctm.append(ifctm1)
        lftfm.append(iftfm1)
    a1=lastwo_central_moment(lctd)
    a2=lastwo_central_moment(lcst)
    a3=lastwo_central_moment(lfct)
    a4=lastwo_central_moment(lftf)
    b1=lastwo_central_moment(lctdm)
    b2=lastwo_central_moment(lcstm)
    b3=lastwo_central_moment(lfctm)
    b4=lastwo_central_moment(lftfm)
    standard_USR=a1+a2+a3+a4
    m_weight_USR=b1+b2+b3+b4
    extra_descriptor_USR=standard_USR+m_weight_USR
    return extra_descriptor_USR
#------------------------------------------------------------------------------------------
""" HOW TO USE:
>>> from discriminate.usr import shape_recognition
"""
def shape_recognition(moleculeina,moleculeinb):
    v0=descriptors(moleculeina)
    v1=descriptors(moleculeinb)
    manhattan_distance=sum([np.absolute(a-b) for a,b in zip(v0,v1)])
    sab=1.0/(1.0+manhattan_distance/16.0)
    return sab
#------------------------------------------------------------------------------------------
def make_similarity_matrix(moleculein):
    total_molecules, ndescriptors=len(moleculein), 16
    matrix_of_descriptors=np.zeros(shape=(total_molecules,ndescriptors),dtype=float)
    for index in range(total_molecules):
        matrix_of_descriptors[index][:]=descriptors(moleculein[index])
    matrix_of_similarity=np.zeros(shape=(total_molecules,total_molecules),dtype=float)
    for ii in range(total_molecules):
        matrix_of_similarity[ii][ii]=float(1.0)
        vi=matrix_of_descriptors[ii][:]
        for jj in range(ii+1,total_molecules):
            vj=matrix_of_descriptors[jj][:]
            lista=[np.absolute(a-b) for a,b in zip(vi,vj)]
            sij=1.0/(1.0+sum(lista)/float(ndescriptors))
            matrix_of_similarity[ii][jj]=sij
            matrix_of_similarity[jj][ii]=sij
    return matrix_of_similarity
#------------------------------------------------------------------------------------------
def kick_similar_molecules(moleculein, tol=0.95, silence=0, mylogfile=log_file):
    s=[[imol,xmol.info['e']] for imol,xmol in enumerate(moleculein)]
    t = sorted(s, key=lambda x: float(x[1]))
    indexlist=[x[0] for x in t]
    molx=moleculein.copy()
    moly=[molx[x] for x in indexlist]
    if silence==1:
        matrixs=make_similarity_matrix(moly)
        deij=float(0.0)
    #------------------------------------------
    sim_list, no_sim_list, large_list, egral_list = [], [], [], []
    ##no_sim_list != sim_nosim_svar_list if a\=b\=c
    ##no_sim_list: final list
    ##sim_nosim_svar_list: only for print
    total_molecules=len(moleculein)
    for ii in range(total_molecules):
        if ii not in sim_list: no_sim_list.append(indexlist[ii])
        for jj in range(ii+1,total_molecules):
            namei, namej=moly[ii].info['i'], moly[jj].info['i']
            if silence==0:
                deij=np.absolute(moly[ii].info['e'] - moly[jj].info['e']) #Getting energy
                svar=shape_recognition(moly[ii],moly[jj]) if (deij <= tole) else 0.0
            elif silence==1:
                svar=matrixs[ii][jj]
            if ( svar >= tol ):
                if (jj not in sim_list):
                    sim_list.append(jj)
                    large_list.append([indexlist[jj], indexlist[ii], namej, namei, deij, svar])
            else:
                egral_list.append([indexlist[jj], indexlist[ii], namej, namei, deij, svar])
    #------------------------------------------
    t = sorted(large_list, key=lambda x: int(x[0]))
    if silence==0:
        logfile = open(mylogfile,'a')
        print("\nsimilarity_tol = %4.2f" %(tol), file=logfile)
        for ii, iix in enumerate(t):
            jj=str(ii+1).zfill(5)
            chain="%s %s ... DISCRIMINATED: Sim to %s (delE=%6.4f) (s=%5.3f)"
            tuple=(jj, iix[2], iix[3], iix[4], iix[5])
            print(chain %tuple, file=logfile)
        logfile.close()
    if verbose==1:
        logfile = open(disc_file,'a')
        for iix in t:
            xmol=moleculein[iix[0]].copy()
            print(len(xmol), file=logfile)
            chain="%12.8f %s DISCR1 Sim to %s (delE=%f) (s=%4.2f >= tol=%s)"
            tuple=(xmol.info['e'], iix[2], iix[3], iix[4], iix[5],str(tol))
            print(chain %tuple, file=logfile)
            for iatom in xmol:
                chain="%s %16.9f %16.9f %16.9f"
                ipos = iatom.position
                tuple=(iatom.symbol, ipos[0], ipos[1], ipos[2])
                print(chain %tuple, file=logfile)
            xmol=moleculein[iix[1]].copy()
            print(len(xmol), file=logfile)
            print("%12.8f %s" %(xmol.info['e'], xmol.info['i']), file=logfile)
            for iatom in xmol:
                chain="%s %16.9f %16.9f %16.9f"
                ipos = iatom.position
                tuple=(iatom.symbol, ipos[0], ipos[1], ipos[2])
                print(chain %tuple, file=logfile)
        logfile.close()
        logfile = open('details.txt','a')
        print("\nenergy_tol = %4.2f similarity_tol = %4.2f" %(tole,tol), file=logfile)
        u = sorted(egral_list, key=lambda x: int(x[1]))
        v = sorted(u, key=lambda x: int(x[0]))
        for ii, iix in enumerate(v):
            jj=str(ii+1).zfill(5)
            if iix[4]>tole:
                chain1="%s %s DISIMILAR to %s (delE=%6.4f> %4.2f)"
                tuple1=(jj, iix[2], iix[3], iix[4], tole)
                print(chain1 %tuple1, file=logfile)
            if iix[4]<=tole:
                chain2="%s %s DISIMILAR to %s (delE=%6.4f<=%4.2f)(sij=%6.4f<%4.2f)"
                tuple2=(jj, iix[2], iix[3], iix[4], tole, iix[5], tol)
                print(chain2 %tuple2, file=logfile)
        logfile.close()
    #------------------------------------------
    no_sim_list=sorted(no_sim_list)
    moleculeout=[moleculein[x] for x in no_sim_list]
    return moleculeout
#------------------------------------------------------------------------------------------
def make_similarity_matrix_compare(moleculein, moleculeref):
    ndescriptors=16
    total_molecules1=len(moleculein)
    total_molecules2=len(moleculeref)
    matrix_of_descriptors1=np.zeros(shape=(total_molecules1,ndescriptors),dtype=float)
    matrix_of_descriptors2=np.zeros(shape=(total_molecules2,ndescriptors),dtype=float)
    for imol in range(total_molecules1):
        matrix_of_descriptors1[imol][:]=descriptors(moleculein[imol])
    for imol in range(total_molecules2):
        matrix_of_descriptors2[imol][:]=descriptors(moleculeref[imol])
    matrix_of_similarity=np.zeros(shape=(total_molecules1,total_molecules2),dtype=float)
    for imol in range(total_molecules1):
        vi=matrix_of_descriptors1[imol][:]
        for jmol in range(total_molecules2):
            vj=matrix_of_descriptors2[jmol][:]
            lista=[np.absolute(a-b) for a,b in zip(vi,vj)]
            sij=1.0/(1.0+sum(lista)/float(ndescriptors))
            matrix_of_similarity[imol][jmol]=sij
    return matrix_of_similarity
#------------------------------------------------------------------------------------------
def molin_sim_molref(moleculein, moleculeref, tol=0.95, silence=0):
###HOW TO USE:
###from discriminate.usr import molin_sim_molref
    if silence==1:
        matrixs=make_similarity_matrix_compare(moleculein, moleculeref)
        deij=float(0.0)
    moleculeout, count=[], 0
    for ii,imol in enumerate(moleculein):
        ans=0
        for jj,jmol in enumerate(moleculeref):
            if silence==0: deij=np.absolute(jmol.info['e'] - imol.info['e']) #Getting energy
            if (deij <= tole):
                svar=matrixs[ii][jj] if silence==1 else shape_recognition(imol,jmol)
                if (svar >= tol):
                    ans, count=1, count+1
                    #------------------------------------------
                    if silence==0:
                        logfile = open(log_file,'a')
                        kk=str(count).zfill(5)
                        cadena="%s %s ... DISCRIMINATED: Sim to %s (delE=%6.4f) (s=%5.3f)"
                        tuplex=(kk,imol.info['i'], jmol.info['i'], deij, svar)
                        print(cadena %tuplex, file=logfile)
                        logfile.close()
                    if verbose==1:
                        logfile = open(disc_file,'a')
                        for xmol in [imol, jmol]:
                            print(len(xmol), file=logfile)
                            if xmol.info['i']==imol.info['i']:
                                cadena="%12.8f %s DISCRIMINATED: Sim to %s (delE=%6.4f) (sym=%5.3f >= tol=%4.2f)"
                                tuplex=(xmol.info['e'], imol.info['i'], jmol.info['i'], deij, svar, tol)
                                print(cadena %tuplex, file=logfile)
                            else:
                                print("%12.8f  %s" %(xmol.info['e'], xmol.info['i']), file=logfile)
                            for iatom in xmol:
                                cadena="%s %16.9f %16.9f %16.9f"
                                ipos = iatom.position
                                tuple=(iatom.symbol, ipos[0], ipos[1], ipos[2])                                
                                print(cadena %tuple, file=logfile)
                        logfile.close()
                    #------------------------------------------
                    break
        if ans==0:
            xmol=imol.copy()
            moleculeout.extend([xmol])
    if count==0 and silence==2:
        logfile = open(log_file,'a')
        print("ZERO elements discriminated by Similarity (tol = %1.2f): molin_sim_molref" %(tol), file=logfile)
        logfile.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def compare_with_memory(moleculeinalone, moleculereflist, tol=0.95, silence=0):
###from discriminate.usrp import compare_with_memory
    ans=0
    if moleculereflist==[]: return ans
    else:
        molin=moleculeinalone.copy()
        for jmol in moleculereflist:
            svar=shape_recognition(molin, jmol)
            if ( svar > tol ):
                if silence==0:
                    logfile = open(log_file,'a')
                    print("Similarity found in %s with sym=%5.3f >= tol=%4.2f" %(jmol.info['i'], svar, tol), file=logfile)
                    logfile.close()
                ans=1
                break
        return ans
#------------------------------------------------------------------------------------------

