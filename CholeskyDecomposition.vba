Function Cholesky(x As Object)
'Cholesky decomposition
'original 6/14/2004
'I forgot where I refer from...
    
    'decralations
    Dim i, j, k, II, JJ, N, Im() As Integer
    Dim a(), U() As Double 'U() is upper right triangle matrix
    Dim b1, b2, b4 As Double
    
    N = x.Columns.Count 'Degree of freedom of matrix
    ReDim a(1 To N, 1 To N), U(1 To N, 1 To N)
    
'initialize matrix A() & U()
    For i = 1 To N
        For j = 1 To N
            a(i, j) = x.Cells(i, j).Value
            U(i, j) = 0
        Next j
    Next i

'upper right triangle
    'U11
    b1 = Sqr(a(1, 1))
    U(1, 1) = b1

    'U1j
    For j = 2 To N
        U(1, j) = a(1, j) / b1
    Next j
    
    'Uii & Uij
    For i = 2 To N
        b1 = a(i, i)
        b4 = 0
            'Uii
            For k = 1 To i - 1
                b1 = b1 - U(k, i) ^ 2
            Next k
            b2 = Sqr(b1)
            U(i, i) = b2
            'Uij
            If i < N Then
                For j = i + 1 To N
                    b1 = a(i, j)
                    For k = 1 To i - 1
                        b1 = b1 - U(k, i) * U(k, j)
                    Next k
                    U(i, j) = b1 / b2
                Next j
            End If
    Next i
    
    'Output
    Cholesky = U
    
    End Function
