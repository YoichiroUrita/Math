Function Eigenvec(x, y As Object)
'extraction Eiven vectors by Gauss elimination
'x: The original matrix for Eigen value problem
'y: Eigen values by QR or another methods (2 rows : 1st row is real, 2nd is imaginary)
'I forgot where this is referred to.
'6/12/2014
'
Dim a(), eig(), evector() As Double
Dim i, j, k, kk As Integer
Dim b1, pivotvalue, valueb As Double

N = x.Columns.Count
M = y.Columns.Count
ReDim a(1 To N, 1 To N)
ReDim eig(1 To 2, 1 To 2 * N)
ReDim evector(1 To N, 1 To 2 * N)


For i = 1 To 2
    For j = 1 To N
        eig(i, j) = y(i, j)
    Next j
Next i

For kset = 1 To N
    If icc = 0 Then
        alp = eig(1, kset)
        bet = eig(2, kset)
        If bet = 0 Then
            icc = 0
            For i = 1 To N
                For j = 1 To N
                    a(i, j) = x(i, j)
                Next j
                a(i, i) = a(i, i) - alp
            Next i
            nn = N
        Else
            ReDim a(1 To 2 * N, 1 To 2 * N)
            icc = 1
            For i = 1 To N
                ib = N + i
                For j = 1 To N
                    a(i, j) = x(i, j)
                    a(ib, N + j) = x(i, j)
                Next j
                a(i, i) = a(i, i) - alp
                a(ib, ib) = a(ib, ib) - alp
                a(i, ib) = bet
                a(ib, i) = -bet
            Next i
            nn = 2 * N
        End If
    Else
        icc = 2
    End If
    
    If icc < 2 Then
        evector(kset, 1) = 1#
        If a(1, 1) <> 0 And a(2, 1) = 0 Then
            For j = 1 To nn
                a(2, j) = a(1, j)
            Next j
        End If
        For i = 2 To nn
            evector(kset, i) = -a(i, 1)
        Next i
        For i = 2 To nn
            a(i, 1) = 0
        Next i

        'Gaus
        For k = 2 To nn
            ic = 0
            For i = k To nn
                If a(i, k) <> 0 Then
                    ic = 1
                    kk = i
                    Exit For
                End If
            Next i
    
            If ic = 0 Then
                aa$ = "Unconvergence"
                ix = MsgBox(aa$, 0)
            End If
    
            pivotvalue = a(kk, k)
            For j = k To nn
                b1 = a(kk, j)
                a(kk, j) = a(k, j)
                a(k, j) = b1 / pivotvalue
            Next j
            b1 = evector(kset, kk)
            evector(kset, k) = b1 / pivotvalue

            For i = 2 To nn
                If i <> k Then
                    valueb = a(i, k)
                    For j = k To nn
                        a(i, j) = a(i, j) - valueb * a(k, j)
                    Next j
                    evector(kset, i) = evector(kset, i) - valueb * evector(kset, k)
                End If
            Next i
        Next k
    End If  'Gaus end
    
    If icc = 2 Then
        For i = 1 To nn
            If i < N Then
                evector(kset, i) = evector(kset - 1, i)
            Else
                evector(kset, i) = -evector(kset - 1, i)
            End If
        Next i
        icc = 0
    End If
    
    'normalize max=1
    Max = 0
    For i = 1 To N
        If Abs(Max) < Abs(evector(kset, i)) Then
            Max = evector(kset, i)
        End If
    Next i
    For i = 1 To nn
        evector(kset, i) = evector(kset, i) / Max
    Next i
Next kset

Eigenvec = evector

End Function
