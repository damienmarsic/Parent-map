#!/usr/bin/env python

# parent-map version 1.0
# Author: Damien Marsic, damien.marsic@aliyun.com
# 2020-04-26
# License: GNU General Public v3 (GPLv3)

import argparse
from gooey import Gooey, GooeyParser
import sys
from collections import defaultdict
import webbrowser
from os import system
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
cli=False
if len(list(sys.argv))>1:
    cli=True

def parse_CLI():
    parser=argparse.ArgumentParser(description="Analyze parental contributions to protein sequences")
    parser.add_argument('Input',type=str,help="File containing the sequence(s) to be analyzed")
    parser.add_argument('Parent',type=str,help="File containing the parental sequence(s)")
    parser.add_argument('-o','--Output',type=str,help="Output file prefix")
    parser.add_argument('-n','--NumSeq',type=int,help="Number of sequences to analyze (default: all sequences)")
    parser.add_argument('-m','--MinFragLen',type=int,help="Minimal fragment length (default: 6 for protein, 18 for DNA)")
    parser.add_argument('-v','--MinOverlap',type=int,help="Minimal overlap length (default: 2 for protein, 6 for DNA)")
    parser.add_argument('-s','--MaxNameSize',type=int,default=12,help="Maximal size of sequence names, longer names will be replaced with a number (default: 12)")
    parser.add_argument('-c','--SeqChars',type=int,default=120,help="Number of sequence characters per line (default: 120)")
    parser.add_argument('-l','--LowerCase',default=False,action='store_true',help="Display sequences in lower case (default: upper case)")
    parser.add_argument('-e','--VRSides',type=int,default=1,help="Number of characters to include each side of variable regions (default:1)")
    parser.add_argument('-f','--Overwrite',default=False,action='store_true',help="Force overwrite of existing files (default: exit if file exists)")
    parser.add_argument('-S','--Symbols',type=str,default='. -',help="Symbols for identity, no match, gap (no match is only used if single parent)")
    parser.add_argument('-d','--DisplayResults',default=False,action='store_true',help="Display results in browser (default: yes)")
    return parser.parse_args()

@Gooey(show_restart_button=False,menu=[{'name':'About','items':[{
    'type': 'AboutDialog',
    'menuTitle': 'About Parent map',
    'name': 'Parent map',
    'description': 'Characterize protein sequence variants against possible parental sequences',
    'version': '1.0',
    'copyright': '2020',
    'website': '',
    'developer': 'Damien Marsic',
    'license': 'GNU General Public v3 (GPLv3)'}]}])
def parse_GUI():
    parser=GooeyParser(description="Analyze parental contributions to protein sequences")
    parser.add_argument('Input',widget="FileChooser",type=str,help="File containing the sequence(s) to be analyzed")
    parser.add_argument('Parent',widget="FileChooser",type=str,help="File containing the parental sequence(s)")
    parser.add_argument('-o','--Output',type=str,help="Output file prefix")
    parser.add_argument('-n','--NumSeq',type=int,help="Number of sequences to analyze (default: all sequences)")
    parser.add_argument('-m','--MinFragLen',type=int,help="Minimal fragment length (default: 6 for protein, 18 for DNA)")
    parser.add_argument('-v','--MinOverlap',type=int,help="Minimal overlap length (default: 2 for protein, 6 for DNA)")
    parser.add_argument('-s','--MaxNameSize',type=int,default=12,help="Maximal size of sequence names, longer names will be replaced with a number (default: 12)")
    parser.add_argument('-c','--SeqChars',type=int,default=120,help="Number of sequence characters per line (default: 120)")
    parser.add_argument('-l','--LowerCase',default=False,action='store_true',help="Display sequences in lower case (default: upper case)")
    parser.add_argument('-e','--VRSides',type=int,default=1,help="Number of characters to include each side of variable regions (default:1)")
    parser.add_argument('-f','--Overwrite',default=False,action='store_true',help="Force overwrite of existing files (default: exit if file exists)")
    parser.add_argument('-S','--Symbols',type=str,default='. -',help="Symbols for identity, no match, gap (no match is only used if single parent)")
    parser.add_argument('-d','--DisplayResults',default=True,action='store_false',help="Display results in browser (default: yes)")
    return parser.parse_args()

def check_file(filename,txt,x):
    try:
        open(filename)
    except IOError:
        if x==1:
            print('\n  '+txt+' file "'+filename+'" not found!\n')
        return False
    else:
        if x==2:
            print('\n  File '+filename+' already exists! Use -f option to force overwrite or -o option to enter a new name.\n')
        return True
def check_int(arg,lim,msg):
    if type(arg) is int and arg<lim:
        print('\n  '+msg+'\n')
        return False
    else:
        return True
def match(seq,ref,size):
    global minfrag
    matches=[]
    if not size:
        return matches
    a=0
    b=size
    c=0
    while b<=len(seq):
        while b<=len(seq) and seq[a:b] not in ref[c:]:
            a+=1
            b+=1
        if b>len(seq):
            continue
        while b<=len(seq) and seq[a:b] in ref[c:]:
            b+=1
        b-=1
        d=ref.find(seq[a:b],c)
        if b-a>=minfrag or d-c<=b-a:
            c=d
            matches.append((a,b,c))
            a=b
        else:
            a+=1
        b=a+size
    return matches

def hiscore(str1,str2):
    if len(str1)>len(str2):
        a=str1
        b=str2
    else:
        a=str2
        b=str1
    r=len(a)-len(b)
    q=range(len(a))
    data=list(range(r))
    combine(q,a,b,r,0,data,0)
    return

def combine(q,a,b,r,index,data,i):
    global top,score
    if index==r:
        x=0
        for j in range(r):
            if data[j]==0 or (j>0 and data[j]==data[j-1]+1):
                continue
            if j==0:
                m=0
                n=0
            else:
                m=data[j-1]+1
                n=data[j-1]-j+1
            for k in range(len(a[m:data[j]])):
                if a[m:data[j]][k]==b[n:data[j]-j][k]:
                    x+=1
        for k in range(len(a[data[j]+1:])):
            if a[data[j]+1:][k]==b[data[j]-j:][k]:
                x+=1
        if x>score:
            score=x
            top=tuple(data)
        return
    if i>=len(q):
        return
    data[index]=q[i]
    combine(q,a,b,r,index+1,data,i+1)
    combine(q,a,b,r,index,data,i+1)

def refine(i,a,b,x,y,seq,ref,db,sd):
    global minfrag,top,score,alnfail
    if b-a==y-x:
        p,q='',''
        s=0
        for j in range(b-a):
            if seq[a+j]==ref[x+j]:
                if p=='':
                    p=a+j
                if q!='':
                    db.insert(i+1+s,(q,a+j,seq[q:a+j],x+q-a,ref[x+q-a:x+j]))
                    s+=1
                    q=''
                continue
            if p!='':
                db.insert(i+1+s,(p,a+j,x+p-a))
                s+=1
                p=''
            if q=='':
                q=a+j
        if p!='':
            db.insert(i+1+s,(p,a+j+1,x+p-a))
        elif q!='':
            db.insert(i+1+s,(q,a+j+1,seq[q:a+j+1],x+q-a,ref[x+q-a:x+j+1]))
    elif x==y:
        db.insert(i+1,(a,b,seq[a:b],x))
    elif a==b:
        sd.append((a,y-x,ref[x:y],x))
    else:
        A=[]
        p=min(minfrag,b-a)
        p=max(p,minfrag*2//3)
        while not A and (p>minfrag*2//3 or (abs(b-a-y+x)>3 and ((seqtype=='DNA' and p>6) or (seqtype!='DNA' and p>3)))):
            p-=1
            A=match(seq[a:b],ref[x:y],p)
        if A:
            for k in range(len(A)):
                db.insert(i+1+k,(A[k][0]+a,A[k][1]+a,A[k][2]+x))
        elif not A and max(len(seq[a:b]),len(ref[x:y]))*abs(b-a-y+x)>200:
            alnfail=True
            return db,sd
        else:
            top=''
            score=0
            hiscore(seq[a:b],ref[x:y])
            ins=[]
            for j in range(len(top)):
                if j==0:
                    ins.append((top[0],top[0]))
                    continue
                if top[j]>top[j-1]+1:
                    ins.append((top[j],top[j]))
                else:
                    ins[-1]=(ins[-1][0],top[j])
            if score==0:
                l=max(y-x,b-a)
                m=min(y-x,b-a)
                ins=[(m,l-1)]
            l=0
            if b-a>y-x:
                for j in ins:
                    db.insert(i+1,(a+j[0],a+j[1]+1,seq[a+j[0]:a+j[1]+1],x+j[0]-l))
                    l+=j[1]-j[0]+1
            else:
                for j in ins:
                    sd.append((a+j[0]-l,j[1]-j[0]+1,ref[x+j[0]:x+j[0]+j[1]-j[0]+1],x+j[0]))
                    l+=j[1]-j[0]+1
    return db,sd

if not cli:
    args=parse_GUI()
else:
    args=parse_CLI()
fail=False
top=''
score=0
results=[]

# Check arguments

if not check_file(args.Input,'Sequence',1):
    sys.exit()
if not check_file(args.Parent,'Parent',1):
    sys.exit()
if not check_int(args.NumSeq,1,'There should be at least one sequence to be analyzed!'):
    sys.exit()
if not check_int(args.MinFragLen,1,'Minimal fragment length should be >0!'):
    sys.exit()
if not check_int(args.MinOverlap,0,'Minimal overlap length should be >=0!'):
    sys.exit()
if not check_int(args.MaxNameSize,0,'Maximal size of sequence names can not be less than 0!'):
    sys.exit()
if not check_int(args.SeqChars,10,'There should be at least 10 sequence characters per line!'):
    sys.exit()
if not check_int(args.VRSides,0,'The number of characters each side of variable regions can not be less than 0!'):
    sys.exit()
if args.VRSides>args.SeqChars/3:
    print('\n  The number of characters each side of variable regions can not be greater than a third the umber of characters per line!')
    sys.exit()
if not type(args.Symbols) is str or len(args.Symbols)!=3 or args.Symbols[0]==' ' or len(set(args.Symbols))!=3:
    print('\n  Symbols should be exactly 3 characters, all different, and the first character can not be a blank space!\n')
    sys.exit()
if args.Output:
    x=args.Input[:args.Input.rfind('/')+1]
    y=args.Input[:args.Input.rfind("\\")+1]
    outfile=x+y+args.Output
else:
    x=args.Parent
    x=x[x.rfind('/')+1:x.rfind('.')]
    x=x[x.rfind("\\")+1:]
    outfile=args.Input[:args.Input.rfind('.')]+'-'+x
    if args.NumSeq:
        outfile+='-n'+str(args.NumSeq)
    if args.MinFragLen:
        outfile+='-m'+str(args.MinFragLen)
    if args.MinOverlap:
        outfile+='-v'+str(args.MinOverlap)
    if args.MaxNameSize!=8:
        outfile+='-s'+str(args.MaxNameSize)
    if args.SeqChars!=70:
        outfile+='-c'+str(args.SeqChars)
    if args.LowerCase:
        outfile+='-lc'
    if args.VRSides!=1:
        outfile+='-e'+str(args.VRSides)
if not args.Overwrite and (check_file(outfile+'-stats.txt','',2) or check_file(outfile+'-par.txt','',2) or check_file(outfile+'-aln.txt','',2) or check_file(outfile+'-def.txt','',2)):
    sys.exit()

# Open parental sequence(s)

f=open(args.Parent,'r')
temp=f.read().strip()
f.close()
ref={}
if '>' in temp:
    while temp[0]!='>':
        temp=temp[1:]
    temp=temp.split('>')
else:
    temp=['Parent\n'+temp]
for n in temp:
    y=n.split()
    if not n or n[0] in (' ','\n') or len(y)==1 or (len(y)==2 and y[1].isdigit()):
        continue
    if len(y[0])>args.MaxNameSize:
        print('\n  Parental sequence name too long! Decrease name length or increase max name size for -s option!\n\n')
        fail=True
        break
    if y[0] in ref:
        print('\n  Duplicate name found in parental sequences! All sequences must have a different name!\n\n')
        fail=True
        break
    if y[1].isdigit():
        ref[y[0]]=int(y[1])-1
        x=2
    else:
        ref[y[0]]=0
        x=1
    n=n.split(None,x)[-1]
    n=n.replace(' ','')
    n=n.replace('\n','')
    if args.LowerCase:
        n=n.lower()
    else:
        n=n.upper()
    ref[y[0]]=(ref[y[0]],n)
if fail:
    sys.exit()

# More checks

if len(ref)==0:
    print('\n  No sequence was found in file '+args.Parent+' !\n\n')
    sys.exit()
seqtype=''
for n in ref:
    temp='DNA'
    for x in ref[n][1]:
        if x not in 'ATGCatcg':
            temp='other'
            break
    if seqtype and temp!=seqtype:
        break
    seqtype=temp
if temp!=seqtype:
    print('\n  All parental sequences must be of same type (DNA or protein) !\n\n')
    sys.exit()

# Process input file

f=open(args.Input,'r')
g=open(outfile+'-par.txt','w')
wrote=False
name=''
seq=''
num=0
stats={'Name':[],'Length':[],'Parents':[],'Main':[],'Coverage':[],'Matches':[],'ID%':[],'Identities':[],'Ins_sites':[],'Ins':[],'Del_sites':[],'Dels':[],'Subs':[],'Other':[]}
ID=args.Symbols[0]
NM=args.Symbols[1]
GA=args.Symbols[2]
COMB={}
COMB2={}
while True:
    line=f.readline()
    if name and (not line or line[0]=='>'):
        seq=seq.replace(' ','')
        num+=1
        seqdel=[]
        seqdel2=[]
        if args.LowerCase:
            seq=seq.lower()
        else:
            seq=seq.upper()
        if num==1:
            temp='DNA'
            for x in seq:
                if x not in 'ATGCatcg':
                    temp='other'
                    break
            if temp!=seqtype:
                print('\n  Sequences to be analyzed must be of same type (DNA or protein) as parental sequences !\n\n')
                fail=True
                break
            if args.MinFragLen:
                minfrag=args.MinFragLen
            else:
                if seqtype=='DNA':
                    minfrag=18
                else:
                    minfrag=6
            if args.MinOverlap:
                minov=args.MinOverlap
            else:
                if seqtype=='DNA':
                    minov=6
                else:
                    minov=2

# Match sequence fragments to parental sequences

        comp=defaultdict(list)
        for m in ref:
            n=ref[m][1]
            comp[m].extend(match(seq,n,minfrag))

# Remove irrelevant matches

        temp=set()
        for m in comp:
            for n in comp[m]:
                temp.add(n)
        a=[]
        while temp:
            a.append(min(temp))
            temp.remove(min(temp))
            x=[]
            for m in comp:
                for n in comp[m]:
                    if n[:2]==a[-1][:2]:
                        x.append(m)
            if len(x)>1:
                y=0
                z=''
                for n in range(len(x)):
                    b=0
                    for i in comp[x[n]]:
                        b+=i[1]-i[0]
                    if b<y:
                        z=x[n]
                    elif b>y and y>0:
                        z=x[n-1]
                        y=b
                    else:
                        y=b
                    if z:
                        for j in comp[z]:
                            if j[:2]==a[-1][:2]:
                                comp[z].remove(j)
                                break
            if len(a)>=2:
                x=''
                if a[-2][1]>a[-1][1] or (a[-2][1]==a[-1][1] and a[-2][0]<a[-1][0]):
                    x=a[-1]
                elif a[-2][0]==a[-1][0] and a[-1][1]>a[-2][1]:
                    x=a[-2]
                if x:
                    for m in comp:
                        for n in comp[m]:
                            if n[:2]==x[:2]:
                                comp[m].remove(n)
                                break
                    if x==a[-2]:
                        del a[-2]
                    else:
                        del a[-1]
            if len(a)<3:
                continue
            if a[0][1]-a[2][0]>=minov and (a[2][0]>a[0][0] or a[2][1]>a[0][1]):
                for m in comp:
                    if a[1] in comp[m]:
                        comp[m].remove(a[1])
                del a[1]
            else:
                del a[0]
        temp=set()
        for m in comp:
            for n in comp[m]:
                temp.add(n)
        a=[]
        while temp:
            a.append(min(temp))
            temp.remove(min(temp))
            if len(a)<3:
                continue
            if (a[1][0]>a[0][1] and a[1][1]<a[2][1]) or (a[1][0]>a[0][0] and a[1][1]<a[2][0]):
                for m in comp:
                    if a[0] in comp[m] and a[2] in comp[m]:
                        for n in comp:
                            if n==m:
                                continue
                            if a[1] in comp[n]:
                                comp[n].remove(a[1])
                        del a[1]
                        break
            del a[0]
        for m in comp:
            if len(comp[m])==1:
                x=comp[m][0]
                for n in comp:
                    if n==m:
                        continue
                    for i in range(len(comp[n])):
                        if comp[n][i][1]>=x[0] and len(comp[n])>i+1 and x[1]>=comp[n][i+1][0] and comp[n][i+1][0]==comp[n][i][1]+1:
                            comp[m].remove(x)
                            break
        x=list(comp.keys())
        for m in x:
            if not comp[m]:
                del comp[m]

# Refine coverage of unmatched regions flanked by same parent or end regions
# Also save relevant info for stats
        temp=set()
        for m in comp:
            for n in comp[m]:
                temp.add(n)
        temp=sorted(temp)
        i=-1
        while i<len(temp) and len(temp)>0:
            if i>=0 and i+1<len(temp):
                for n in comp:
                    if temp[i] in comp[n] and temp[i+1] in comp[n]:
                        break
                else:
                    i+=1
                    continue
            if i==-1 and temp[0][0]>0:
                a=0
                x=0
                for n in comp:
                    if temp[0] in comp[n]:
                        break
            elif i==-1 and temp[0][0]==0:
                i+=1
                continue
            else:
                a=temp[i][1]
                if len(temp[i])==3:
                   x=temp[i][2]+a-temp[i][0]
                elif len(temp[i])==5:
                   x=temp[i][3]+a-temp[i][0]
                elif len(temp[i])==4:
                   x=temp[i][3]
            if i+1==len(temp) and temp[i][1]<len(seq):
                for n in comp:
                    if temp[i] in comp[n]:
                        break
                b=len(seq)
                y=len(ref[n][1])
            elif i+1==len(temp) and temp[i][1]==len(seq):
                i+=1
                continue
            else:
                b=temp[i+1][0]
                if len(temp[i+1])==3:
                    y=temp[i+1][2]
                else:
                    y=temp[i+1][3]
            if a==b and (x==y or (i+1<len(temp) and len(temp[i])!=len(temp[i+1])) or a in [k for k,_,_,_ in seqdel]):
                i+=1
                continue
            for k,l,_,m in seqdel:
                if a==k:
                    x+=l
                if a<k<=b:
                    b=k
                    y=m
                    break
            m=0
            while m<len(temp):
                if temp[m] not in comp[n]:
                    temp.remove(temp[m])
                else:
                    m+=1
            alnfail=False
            temp,seqdel=refine(i,a,b,x,y,seq,ref[n][1],temp,seqdel)
            for m in temp:
                if not m in comp[n]:
                    comp[n].append(m)
            comp[n]=sorted(comp[n])
            temp=set()
            for m in comp:
                for k in comp[m]:
                    temp.add(k)
            temp=sorted(temp)
            if alnfail:
                print('\n  Not enough sequence similarity between '+name+' and '+n+'!\n  Aborting refinement...')
                break

# Write parental maps to file

        g.write('Parental composition of sequence '+name)
        if len(name)>args.MaxNameSize:
            print("\n  Sequence name '"+name+"' replaced with ",end='')
            name='S'+str(num)
            if len(name)>args.MaxNameSize:
                print('\n\n  Maximal name size is too short ! Choose a larger value with option -s\n\n')
                fail=True
                break
            g.write('\nSequence name too long, replaced with '+name+' (the number indicates sequence position in the input file)')
            print("'"+name+"'")
        z=0
        q=0
        MNS=args.MaxNameSize
        b=args.SeqChars
        B=0
        C=0
        while z<len(seq):
            if q+10<=len(seq):
                g.write('\n\n')
                q+=10
                if seqdel:
                    i=0
                    for n in seqdel:
                        if n[0]<q-1 and n[0]+i>q-11:
                            i+=n[1]
                            j=n[0]
                        if i>0 and j!=n[0]:
                            break
                    if i>0 and (B!=0 or j>=z):
                        g.write(' '*(i-B))
                if MNS+1-z+q-len(str(q))>MNS:
                    g.write(' '*(MNS+1-z+q-len(str(q)))+str(q))
                else:
                    g.write(' '*(MNS+1-z+q))
            else:
                g.write('\n')
            i=0
            while q+10+i<=z+b and q+10<=len(seq):
                q+=10
                if seqdel:
                    j=0
                    for n in seqdel:
                        if n[0]<q-1 and n[0]+j>q-11:
                            j+=n[1]
                    if j>0:
                        i+=j
                        if q+i<=z+b:
                            g.write(' '*j)
                if q+i<=z+b:
                    g.write('{:>10}'.format(str(q)))
                else:
                    q-=10
            g.write('\n'+name+' '*(MNS+1-len(name)))
            wrote=True
            A=0   # cumulated number of sequence deletions in current line
            y=0   # seqdel index
            C=B
            while seqdel and y<len(seqdel) and z>seqdel[y][0]:
                y+=1
            for n in range(z,z+b):
                if n-A>=len(seq):
                    break
                if seqdel and y<len(seqdel) and C==seqdel[y][1]:
                    y+=1
                    C=0
                if not seqdel:
                    g.write(seq[n])
                if seqdel and y<len(seqdel) and C<seqdel[y][1] and n-A>=seqdel[y][0]:
                    g.write(GA)
                    A+=1
                    C+=1
                elif seqdel:
                    g.write(seq[n-A])
            for m in comp:
                n=comp[m]
                x=0
                while x<len(n) and z>=n[x][1]:
                    x+=1
                if x==len(n):
                    continue
                g.write('\n'+m+' '*(MNS+1-len(m)))
                A=0
                y=0
                while seqdel and y<len(seqdel) and z>seqdel[y][0]:
                    y+=1
                for i in range(z,z+b):
                    if i-A==n[x][1]:
                        x+=1
                        if x==len(n):
                            break
                    if seqdel and y<len(seqdel) and B==seqdel[y][1] :
                        y+=1
                        B=0
                    if seqdel and y<len(seqdel) and B<seqdel[y][1] and i-A>=seqdel[y][0]:
                        g.write(seqdel[y][2][B])
                        A+=1
                        B+=1
                    elif i-A<n[x][0]:
                        if len(comp)==1:
                            g.write(NM)
                        else:
                            g.write(' ')
                    elif i-A<n[x][1]:
                        if len(n[x])==4:
                            g.write(GA)
                        elif len(n[x])==5:
                            g.write(n[x][4][i-A-n[x][0]])
                        else:
                            g.write(ID)
            z+=b-A
        g.write('\n\n\n')

# Save to stats and combine data

        stats['Name'].append(name)
        stats['Length'].append(len(seq))
        stats['Parents'].append(len(comp))
        cov=0
        IS=0
        IC=0
        SC=0
        main=''
        for n in comp:
            x=0
            for m in comp[n]:
                if len(m)==3:
                    x+=m[1]-m[0]
                elif len(m)==4:
                    IS+=1
                    IC+=m[1]-m[0]
                elif len(m)==5:
                    SC+=m[1]-m[0]
            if x>cov:
                cov=x
                main=n
        stats['Main'].append(main)
        stats['Coverage'].append(round(cov/len(seq)*100,2))
        stats['Matches'].append(cov)
        stats['Ins_sites'].append(IS)
        stats['Ins'].append(IC)
        DS=0
        DC=0
        for n in seqdel:
            DS+=1
            DC+=n[1]
        stats['Del_sites'].append(DS)
        stats['Dels'].append(DC)
        stats['Subs'].append(SC)
        OU=0
        y=0
        temp=set()
        for m in comp:
            for n in comp[m]:
                temp.add(n)
        temp=sorted(temp)
        for n in temp:
            if n[0]>y:
                OU+=n[0]-y
                comp[0].append((y,seq[y:n[0]]))
            y=n[1]
        OU+=len(seq)-y
        if len(seq)!=y:
            comp[0].append((y,seq[y:]))
        stats['Other'].append(OU)
        for n in seqdel:
            comp[-1].append(n)
        COMB[name]=comp

# map against main parent for variants from multiple parents

        k=len(comp)
        if 0 in comp:
            k-=1
        if -1 in comp:
            k-=1
        if k>1:
            temp=match(seq,ref[main][1],minfrag)
            i=-1
            while i<len(temp):
                temp=sorted(temp)
                if i==-1:
                    a=0
                    x=0
                else:
                    a=temp[i][1]
                    if len(temp[i])==3:
                       x=temp[i][2]+a-temp[i][0]
                    elif len(temp[i])==5:
                       x=temp[i][3]+a-temp[i][0]
                    elif len(temp[i])==4:
                       x=temp[i][3]
                if i+1==len(temp):
                    b=len(seq)
                    y=len(ref[main][1])
                else:
                    b=temp[i+1][0]
                    if len(temp[i+1])==3:
                        y=temp[i+1][2]
                    else:
                        y=temp[i+1][3]
                if a==b and (x==y or (i+1<len(temp) and len(temp[i])!=len(temp[i+1])) or a in [k for k,_,_,_ in seqdel2]):
                    i+=1
                    continue
                for k,l,_,m in seqdel2:
                    if a==k:
                        x+=l
                    if a<k<=b:
                        b=k
                        y=m
                        break
                alnfail=False
                temp,seqdel2=refine(i,a,b,x,y,seq,ref[main][1],temp,seqdel2)
                if alnfail:
                    print('\n  Not enough sequence similarity between '+name+' and '+main+'!\n  Aborting alignment against main parent...')
                    id=0
                    break
            else:
                for n in seqdel2:
                    temp.append(n)
                temp=sorted(temp)
                COMB2[name]={main:temp}
                id=0
                for m in temp:
                    if len(m)==3:
                        id+=m[1]-m[0]
        else:
            id=cov
        stats['ID%'].append(round(id/len(seq)*100,2))
        stats['Identities'].append(id)

# Next input sequence

    if not line:
        break
    if not line.strip():
        continue
    if line[0]=='>':
        if num==args.NumSeq:
            break
        name=line[1:].split()[0]
        seq=''
        continue
    seq+=line.strip()
    if not name:
        name="Seq"
g.close()
if wrote:
    print('\n  Parental maps saved into file: '+outfile+'-par.txt\n')
    results.append('par')
f.close()
if fail:
    sys.exit()

# Statistics

stats=pd.DataFrame(stats)
stats=stats.set_index('Name')
stats=stats.sort_values(by=['Parents','Main','Ins_sites'])
if wrote:
    g=open(outfile+'-stats.txt','w')
    g.write(str(stats)+'\n\n')
    g.close()
    print('  Stats saved into file: '+outfile+'-stats.txt\n')
    results.append('stats')

# Sequence definitions

g=open(outfile+'-def.txt','w')
wrote=False
for n in COMB:
    A=0
    temp=[]
    for m in COMB[n]:
        temp.extend(COMB[n][m])
    temp=sorted(temp)
    for m in temp:
        for k in m:
            if not str(k).isdigit():
                if len(k)>A:
                    A=len(k)
    g.write('Variant name: '+n+'\n')
    g.write('Variant region  Parent/feature  Parent region  Variant sequence'+' '*(max(A,16)-16)+'  Parent sequence\n')
    for m in temp:
        a=[k for k in COMB[n] if m in COMB[n][k]][0]
        x=str(m[0]+1)
        z=''
        if a==-1:
            x=str(m[0])
            y=str(m[0]+1)
            z='deletion'
            w=str(m[3]+X)
            if len(m[2])!=1:
                w+='-'+str(m[3]+X+len(m[2])-1)
            w=' '*(15-len(w))+w+' '*(4+max(A,16))+m[2]
        elif a==0:
            y=str(m[0]+1+len(m[1])-1)
            z='unmatched'
            w=' '*17+m[1]
        elif a!=-1 and a!=0:
            X=ref[a][0]
            if X<1:
                X=1
            y=str(m[1])
        if len(m)==3:
            z=a
            w=str(m[2]+X)
            if m[1]!=m[0]+1:
                w+='-'+str(m[2]+X+m[1]-m[0]-1)
            w=' '*(15-len(w))+w
        elif len(m)==5:
            z='substitution'
            w=str(m[3]+X)
            if m[1]!=m[0]+1:
                w+='-'+str(m[3]+X+m[1]-m[0]-1)
            w=' '*(15-len(w))+w+'  '+m[2]+' '*(max(A,16)+2-len(m[2]))+m[4]
        elif len(m)==4 and a!=-1:
            z='insertion'
            w=str(m[3]+X-1)+'/'+str(m[3]+X)
            w=' '*(15-len(w))+w+'  '+m[2]
        if x==y:
            q=x
        else:
            q=x+'-'+y
            if a==-1:
                q=x+'/'+y
        g.write(' '*(14-len(q))+q+' '*(16-len(z))+z+w+'\n')
        wrote=True
    g.write('\n')
g.close()
if wrote:
    print('  Sequence definitions saved into file: '+outfile+'-def.txt\n')
    results.append('def')

# Alignments

VRS=args.VRSides
g=open(outfile+'-aln.txt','w')
wrote=False
for u in range(2):
    if (u==1 and len(list(stats.loc[stats['Parents']>1].index))==0) or (u==0 and len(list(stats.loc[stats['Parents']==1].index))==0):
        continue
    g.write('\nAlignment of variant sequences against their main parental sequence')
    if len(ref)>1:
        if u==0:
            g.write(' - single parent\n')
        else:
            g.write(' - multiple parents\n')
    for n in ref:
        temp={}
        if u==0:
            m=list(stats.loc[(stats['Parents']==1) & (stats['Main']==n)].index)
        else:
            m=list(stats.loc[(stats['Parents']>1) & (stats['Main']==n)].index)
        if not m:
            continue
        if u==0:
            for x in m:
                y=set()
                for i in COMB[x]:
                    for j in COMB[x][i]:
                        y.add(j)
                temp[x]=sorted(y)
        else:
            for x in m:
                if x not in COMB2:
                    continue
                y=set()
                for i in COMB2[x]:
                    for j in COMB2[x][i]:
                        y.add(j)
                temp[x]=sorted(y)
        X=ref[n][0]
        if X<1:
            X=1
        A=[]    # A: list of regions to display, from variations in each seq / format: ref coords, max size
        for x in temp:  # x: sequence name
            y=0         # y: list index
            if len(temp[x][y])==2:
                y+=1
            if len(temp[x][y])==3 and temp[x][y][2]!=0:
                A.append((0,temp[x][y][2],max(temp[x][y][2],temp[x][y][0])))
            elif len(temp[x][y])!=3 and temp[x][y][3]!=0:
                A.append((0,temp[x][y][3],max(temp[x][y][3],temp[x][y][0])))
            elif len(temp[x][y])==4 and temp[x][y][3]==0:
                A.append((0,0,len(temp[x][y][2])))
            while y<len(temp[x])-1:
                y+=1
                if len(temp[x][y])==2:  # Unmatched region
                    if y+1==len(temp[x]):
                        z=len(ref[n][1])
                    else:
                        z=temp[x][y+1][2]
                    A.append((temp[x][y-1][2]+temp[x][y-1][1]-temp[x][y-1][0],z,z-temp[x][y-1][2]+temp[x][y-1][1]-temp[x][y-1][0]))
                elif len(temp[x][y])==5:                                                         # Substitution
                    A.append((temp[x][y][3],temp[x][y][3]+len(temp[x][y][4]),len(temp[x][y][4])))
                elif len(temp[x][y])==4 and len(temp[x][y][2])==temp[x][y][1]-temp[x][y][0]:     # Insertion
                    A.append((temp[x][y][3],temp[x][y][3],len(temp[x][y][2])))
                elif len(temp[x][y])==4 and len(temp[x][y][2])==temp[x][y][1]:                   # Deletion
                    A.append((temp[x][y][3],temp[x][y][3]+temp[x][y][1],temp[x][y][1]))
            z=0
            if A and A[-1][1]!=len(ref[n][1]) and len(temp[x][y])==5:
                z=temp[x][y][3]+temp[x][y][1]-temp[x][y][0]
            elif len(temp[x][y])==3 and temp[x][y][2]+temp[x][y][1]-temp[x][y][0]!=len(ref[n][1]):
                z=temp[x][y][2]+temp[x][y][1]-temp[x][y][0]
            if z:
                A.append((z,len(ref[n][1]),len(ref[n][1])-z))
        A=sorted(A)
        C={}                       # list of insertions
        for i in range(len(A)):
            x=max(0,A[i][0]-VRS)
            y=min(len(ref[n][1]),A[i][1]+VRS)
            z=A[i][2]+A[i][0]-x+y-A[i][1]
            if A[i][0]==A[i][1]:
                if A[i][0] in C:
                    C[A[i][0]]=max(A[i][2],C[A[i][0]])
                else:
                    C[A[i][0]]=A[i][2]
            A[i]=(x,y,z)
        if A:
            B=[A[0]]
        else:
            B=[]
        for x in A[1:]:
            y=B[-1]
            if x[1]>=y[0] and y[1]>=x[0]:
                B.remove(y)
                z=0
                for i in C:
                    if y[0]<=i<=x[1]:
                        z+=C[i]
                if z:
                    z+=x[1]-y[0]
                else:
                    z=x[2]+y[2]-y[1]+x[0]
                B.append((y[0],x[1],z))
            else:
                B.append(x)

    # Write alignment to file

        i=0
        j=0
        if len(B)==1:
            j=-1
        while j<len(B)-1:
            L=0
            if j!=0:
                i=j+1
            if len(B)==1:
                i=0
            while j<len(B)-1:
                x=max(len(str(B[j+1][0]+X)),B[j+1][2])+1
                if L+x>args.SeqChars and len(B)>1:
                    break
                L+=x
                j+=1
            if i==j+1 or x>args.SeqChars:
                print('\n  Failed alignment for parent: '+n+'! Try with different values for -c and -e options...\n')
                break
            g.write('\n\n'+' '*MNS)
            wrote=True
            for x in range(i,j+1):
                g.write(' '+str(B[x][0]+X)+' '*(max(0,B[x][2]-len(str(B[x][0]+X)))))
            g.write('\n'+n+' '*(1+MNS-len(n)))
            for x in range(i,j+1):
                w=''
                c=B[x][0]
                for p in sorted(C):
                    if B[x][0]<=p<=B[x][1]:
                        w+=ref[n][1][c:p]+GA*C[p]
                        c=p
                w+=ref[n][1][c:B[x][1]]
                g.write(w+' '*(1+max(0,len(str(B[x][0]+X))-len(w))))
            for m in temp:
                g.write('\n'+m+' '*(1+MNS-len(m)))
                for x in range(i,j+1):
                    z=0
                    w=''
                    k=B[x][0]
                    r=1
                    for y in temp[m]:
                        if len(y)==3:
                            c=y[2]
                            z=y[2]+y[1]-y[0]
                        elif len(y)>3:
                            c=y[3]
                        else:
                            c=z
                        if not k<=c<=B[x][1]:
                            continue
                        for p in sorted(C):
                            if k<p<c or (p==c and len(y)!=4) or (p==c and len(y)==4 and len(y[2])==y[1] and y[0]!=0) or (p==0 and p<c and y[0]==0):
                                w+=ID*(p-k)+GA*C[p]*r
                                r=1
                                k=p
                            if p>c:
                                break
                        if len(y)==2:
                            w+=ID*(c-k)+y[1]
                            k=len(y[1])+c
                        elif len(y)==3:
                            w+=GA*(y[2]-k)
                            k=y[2]
                        elif len(y)==5:
                            w+=ID*(c-k)
                            z=0
                            for s in sorted(C):
                                if c<s<c+len(y[2]):
                                    w+=y[2][z:s-c]+GA*C[s]
                                    z=s-c
                            if z:
                                w+=y[2][z:]
                            else:
                                w+=y[2]
                            k=len(y[2])+c
                        elif len(y)==4 and len(y[2])==y[1]-y[0]:
                            w+=ID*(c-k)+y[2]
                            k=c
                            r=0
                            if c in C:
                                 w+=GA*(C[c]-len(y[2]))
                        elif len(y)==4 and len(y[2])==y[1]:
                            w+=ID*(c-k)+GA*y[1]
                            k=c+y[1]
                    for p in sorted(C):
                        if k<p<B[x][1]:
                            w+=ID*(p-k)+GA*C[p]
                            k=p
                    if x<len(B)-1 or B[x][1]!=len(ref[n][1]):
                        w+=ID*(B[x][1]-k)
                    elif x==len(B)-1 and B[x][1]==len(ref[n][1]):
                        if len(y)==3:
                            z=y[2]+y[1]-y[0]-k
                        else:
                            z=0
                        q=B[x][1]-k-z
                        w+=ID*z+GA*q
                    if x==len(B)-1 and B[x][1] in C and y[-1]!=B[x][1]:
                        w+=GA*C[B[x][1]]
                    g.write(w+' '*(1+max(0,len(str(B[x][0]+X))-len(w))))
        g.write('\n')
    g.write('\n\n')
g.write('\n')
g.close()
if wrote:
    print('  Alignments saved into file: '+outfile+'-aln.txt\n')
    results.append('aln')
if not args.DisplayResults:
    sys.exit()
for x in ("chrome","firefox","iexplore","opera"):
    for n in results:
        URL='file:///'+outfile+'-'+n+'.txt'
        try:
            system("start "+x+" "+URL)
        except:
            break
    else:
        break
else:
    for x in ("chrome","firefox","iexplore","opera"):
        for n in results:
            URL='file:///'+outfile+'-'+n+'.txt'
            try:
                webbrowser.get(x).open(URL)
            except:
                break
        else:
            break
    else:
        for n in results:
            URL='file:///'+outfile+'-'+n+'.txt'
            webbrowser.open_new_tab(URL)
