# Code-Templates

## DSU

  ```
 using namespace std;
       class DisjointSet {
          vector<int> rank, parent, size;
      public:
          DisjointSet(int n) {
              rank.resize(n + 1, 0);
              parent.resize(n + 1);
              size.resize(n + 1);
              for (int i = 0; i <= n; i++) {
              parent[i] = i;
              size[i] = 1;
          }
      }
  
      int findUPar(int node) {
          if (node == parent[node])
              return node;
          return parent[node] = findUPar(parent[node]);
      }
      void unionBySize(int u, int v) {
          int ulp_u = findUPar(u);
          int ulp_v = findUPar(v);
          if (ulp_u == ulp_v) return;
          if (size[ulp_u] < size[ulp_v]) {
              parent[ulp_u] = ulp_v;
              size[ulp_v] += size[ulp_u];
          }
          else {
              parent[ulp_v] = ulp_u;
              size[ulp_u] += size[ulp_v];
          }
      }
  };
```

## Upper bound code
```
int upperBound(vector<int> &arr, int x, int n) {
    int low = 0, high = n - 1;
    int ans = n;

    while (low <= high) {
        int mid = (low + high) / 2;
        // maybe an answer
        if (arr[mid] > x) {
            ans = mid;
            //look for smaller index on the left
            high = mid - 1;
        }
        else {
            low = mid + 1; // look on the right
        }
    }
    return ans;
}

```

## Lower Bound
```
int lowerBound(vector<int> arr, int n, int x) {
    int low = 0, high = n - 1;
    int ans = n;

    while (low <= high) {
        int mid = (low + high) / 2;
        // maybe an answer
        if (arr[mid] >= x) {
            ans = mid;
            //look for smaller index on the left
            high = mid - 1;
        }
        else {
            low = mid + 1; // look on the right
        }
    }
    return ans;
}
```
## Seive - find all primes till N
```
vector<int> findAllPrimes(int n) {
    vector<int> prime(n + 1, 1);
    prime[0] = prime[1] = 0; 
    for (int i = 2; i <= sqrt(n); ++i) {
        if (prime[i] == 1) {
            for (int j = i * i; j <= n; j += i) {
                // Mark multiples of prime
                // numbers as not prime
                prime[j] = 0; 
            }
        }
    }
    
    vector<int> ans;
    for (int i = 2; i <= n; ++i) {
        if (prime[i] == 1) {
            ans.push_back(i);
        }
    }
    return ans;
}
```
## check prime

```
From above code
```



