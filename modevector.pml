run modevectors.py
modevector sys1, sys2
create prot1, sys1 and not resn MLT+LIG+NDP
create prot2, sys2 and not resn MLT+LIG+NDP
color raspberry, prot1
color grey80, prot2
show cartoon, prot1
hide lines, prot1
show cartoon, prot2
hide lines, prot2

create prot1_chaA, prot1 and chain A
create prot1_chaB, prot1 and chain B
create prot2_chaA, prot2 and chain A
create prot2_chaB, prot2 and chain B

select E1, resi 7-10
select E2, resi 28-31
select E3, resi 49-52
select E4, resi 71-74
select E5, resi 93-96
select H1, resi 16-25
select H2, resi 37-46
select H3, resi 60-65
select H4, resi 84-90

select E6, resi 146-150
select E7, resi 169-173
select E8, resi 201-204
select E9, resi 228-232
select E10, resi 253-258
select E11, resi 277-279

select btmBeta, E1 or E2 or E3 or E4 or E5
select topBeta, E6 or E7 or E8 or E9 or E10 or E11

select H5, resi 101-118
select H6, resi 120-128
select H7, resi 154-164
select H8, resi 179-183
select H9, resi 192-194
select H10, resi 218-223
select H11, resi 241-250

select H12, resi 288-305

select L1, resi 53-59
select L2, resi 75-82
select L3, resi 174-178
select L4, resi 205-217
select L5, resi 233-240
select L6, resi 259-276
select L7, resi 130-139
