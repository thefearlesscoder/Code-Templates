# Code-Templates

## DSU

  ```
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
bool checkprime(int n)
{
    if (n <= 1)
        return false;
    for (int i = 2; i * i <= n; i++)
    {
        if (n % i == 0)
            return false;
    }
    return true;
}
```

## Modular exponentiation

```
int exp(int x, int n, int m) {
	assert(n >= 0);
	x %= m;  // note: m * m must be less than 2^63 to avoid ll overflow
	int res = 1;
	while (n > 0) {
		if (n % 2 == 1) { res = res * x % m; }
		x = x * x % m;
		n /= 2;
	}
	return res;
}
```

## insert elements in a Vector, during recursive call

```
vector<int> solve(int idx, int prev,vector<int>& nums, vector<vector<vector<int>>> &dp){
        if(idx >= n)return {};
        if(dp[idx][prev+1].size())return dp[idx][prev+1];

        vector<int> take, notTake;
        if(prev==-1 || nums[idx]%nums[prev] == 0){
            take = solve(idx+1, idx,nums, dp);
            take.insert(take.begin(), nums[idx]);
        }
        notTake = solve(idx+1, prev,nums, dp);

        if(take.size() > notTake.size()){
            return dp[idx][prev+1] = take;
        }else{
            return dp[idx][prev+1] = notTake;
        }
    }

```

## Custom sorting in a vector
```
	sort(edgeList.begin(), edgeList.end(), [](const vector<int>& a, const vector<int>& b) {
            return a[2] < b[2];
        });
```

## Convert Character to String 
```
	char ch = 'a';
	string temp = string(1, ch);
```

## Next Greater Element
```
vector<int> nextGreaterElement(const vector<int>& nums) {
    int n = nums.size();
    vector<int> nge(n, -1);
    stack<int> st;
    for (int i = 2 * n - 1; i >= 0; i--) {
        while (!st.empty() && st.top() <= nums[i % n]) {
            st.pop();
        }
        if (i < n && !st.empty()) {
            nge[i] = st.top();
        }
        st.push(nums[i % n]);
    }
    return nge;
}
```

## Next Smaller Element
```
vector<int> nextSmallerElement(const vector<int>& nums) {
    int n = nums.size();
    vector<int> nse(n, -1);
    stack<int> st;
    for (int i = n - 1; i >= 0; i--) {
        while (!st.empty() && st.top() >= nums[i]) {
            st.pop();
        }
        if (!st.empty()) {
            nse[i] = st.top();
        }
        st.push(nums[i]);
    }
    return nse;
}
```

## Previous Greater Element
```
vector<int> previousGreaterElement(const vector<int>& nums) {
    int n = nums.size();
    vector<int> pge(n, -1);
    stack<int> st;
    for (int i = 0; i < n; i++) {
        while (!st.empty() && st.top() <= nums[i]) {
            st.pop();
        }
        if (!st.empty()) {
            pge[i] = st.top();
        }
        st.push(nums[i]);
    }
    return pge;
}
```

## Previous Smaller Element
```
vector<int> previousSmallerElement(const vector<int>& nums) {
    int n = nums.size();
    vector<int> pse(n, -1);
    stack<int> st;
    for (int i = 0; i < n; i++) {
        while (!st.empty() && st.top() >= nums[i]) {
            st.pop();
        }
        if (!st.empty()) {
            pse[i] = st.top();
        }
        st.push(nums[i]);
    }
    return pse;
}
```
