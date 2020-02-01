#! /usr/bin/env python
# -*- coding: utf-8 -*-

with open('what.list') as ff:
    what_list = ff.read().splitlines()

objs = {}
objcnt = {}
for what_line in what_list:
    if len(what_line.split()) < 4:
        continue
    fname = what_line.split()[0]
    obj = ' '.join(what_line.split()[3:])
    name = what_line.split()[3]
    seg = what_line.split()[4]
    if '[A]' in obj or name not in objcnt:
        cnt = objcnt.get(name, 0) + 1
        vals = objs.get(name, {})
        objcnt[name] = cnt
        vals[cnt] = [fname]
        objs[name] = vals
    else:
        cnt = objcnt[name]
        objs[name][cnt].append(fname)

outlist = []
arepairs = False
for k, v in objs.items():
    if "STD-" not in k:
        for kk, w in v.items():
            objstr = "%s" % k
            if len(w) == 2:
                arepairs = True
                objstr += " %s %s" % (w[0], w[1])
                outlist.append(objstr)

if arepairs:
    outlist.sort()
    f = open("abpairs.tab", "w")
    for outstr in outlist:
        f.write(outstr+'\n')
    f.close()
