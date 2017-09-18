#define RED(a,i,j)  (unsigned int)a[i*width*3 + j*3 + 0]
#define GREEN(a,i,j)  (unsigned int)a[i*width*3 + j*3 + 1]
#define BLUE(a,i,j)  (unsigned int)a[i*width*3 + j*3 + 2]
#define ILUM(a,i,j) (double)(RED(a,i,j) + GREEN(a,i,j) + BLUE(a,i,j)) / 3
#define DEBUG(x) cerr << #x << " = " << x << endl;
#define DEBUGM(a,n,m)  cout << #a << " = " << endl;\
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < m; j++) {\
            cout << a[i][j] << " ";\
        }\
        cout << endl; \
    } \
    cout << endl;\

#define DEBUGV(a,n)  cout << #a << " = " << endl;\
    for(int i = 0; i < n; i++) {\
        cout << a[i] << " ";\
    }\
    cout << endl;
