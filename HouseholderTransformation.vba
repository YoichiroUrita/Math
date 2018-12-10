Function Householder(x As Object)

'Hessenberg Matrix by Householder Transformation
'You can use this function as Excel cell function
'6/6/2014

'
Dim i, j, k, M, II, JJ, N As Integer
Dim a(), buf(), uvec(), qmat() As Double
Dim b1, b2, b3, sigma, c1 As Double

N = x.Columns.Count
ReDim a(1 To N, 1 To N)
ReDim buf(1 To N, 1 To N)
ReDim uvec(1 To N)
ReDim qmat(1 To N, 1 To N)

'initialize
    For i = 1 To N
        For j = 1 To N
            a(i, j) = x.Cells(i, j).Value
            buf(i, j) = 0
            qmat(i, j) = 0
        Next j
    Next i
    
For k = 1 To N - 2
    'U-vector
    b1 = 0
    For i = k + 1 To N
        b1 = b1 + a(i, k) ^ 2
    Next i
    b2 = Sqr(b1)
    b3 = a(k + 1, k)
    If b3 >= 0 Then
        sigma = b2
    Else
        sigma = -b2
    End If
    
    uvec(k + 1) = b3 + sigma
    c1 = sigma * uvec(k + 1)
    uvec(k) = 0
    
    For i = k + 2 To N
        uvec(i) = a(i, k)
    Next i
    
    'Q-matrix

    'initialize
    For i = 1 To N
        For j = 1 To N
            qmat(i, j) = 0
        Next j
    Next i
    For i = 1 To N
        qmat(i, i) = 1
    Next i
    For i = k To N
        For j = k To N
            qmat(i, j) = qmat(i, j) - uvec(i) * uvec(j) / c1
        Next j
    Next i
    
    'calculate QA   Q-Matrix is symmetric and same as inverse matrix
    For i = 1 To k
        For j = k To N
            buf(i, j) = a(i, j)
        Next j
    Next i
    
    For i = k + 1 To N
        For j = k To N
            b1 = 0
            For M = k To N
                b1 = b1 + qmat(i, M) * a(M, j)
            Next M
            buf(i, j) = b1
        Next j
    Next i
    
    'multipul QA and Q by right hand
    For j = 1 To k
        For i = k To N
            a(i, j) = buf(i, j)
        Next i
    Next j
    For j = k + 1 To N
        For i = 1 To N
            b1 = 0
            For M = k To N
                b1 = b1 + buf(i, M) * qmat(M, j)
            Next M
            a(i, j) = b1
        Next i
    Next j
Next k
    
'erase computer error
For i = 3 To N
    For j = 1 To i - 2
        a(i, j) = 0
    Next j
Next i

Householder = a 'if you need Q matrix (qmat) ,replace a to qmat
    
End Function
