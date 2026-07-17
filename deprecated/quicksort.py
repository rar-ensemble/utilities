import numpy as nmp

def qs( dr , inds , left , right ):
  i = left
  j = right
  x = dr[ (left+right) / 2 ]

  while ( ):

    while ( (dr[i] < x ) and ( i < right ) ):
      i += 1
    while ( (x < dr[j])  and ( j > left  ) ):
      j -= 1

    if ( i <= j ):
      y = dr[i]
      dr[i] = dr[j]
      dr[j] = y

      yi = inds[i]
      inds[i] = inds[j]
      inds[j] = yi

      i += 1
      j -= 1
    
    if ( i > j ):
      break

  if ( left < j ):
    qs( dr , inds , left , j )
  if ( i < right):
    qs( dr , inds , i , right)

  return dr , inds
    
