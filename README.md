# Code-Templates - > testting web hook

## check if node u is ancestor of node v
```
// Time of discovery concept using DFS
```

## Fenwik Tree / Binary Indexed Tree
```
```

## Kadane's Algorithm ( maximum sum subarray)
```
int maxSubarraySum(vector<int> &arr) {
    
    // Stores the result (maximum sum found so far)
    int res = arr[0];           
    
    // Maximum sum of subarray ending at current position
    int maxEnding = arr[0];     

    for (int i = 1; i < arr.size(); i++) {
        
        // Either extend the previous subarray or start 
        // new from current element
        maxEnding = max(arr[i], maxEnding + arr[i]);

        // Update result if the new subarray sum is larger
        res = max(res, maxEnding);
    }
    return res;
}
```
## Segment Tree
```
#include<bits/stdc++.h>
using namespace std;

class SegmentTree{
    public:
    vector<int> seg;
    vector<int> lazy;
    vector<int> arr;

    SegmentTree(int n, vector<int> a){
        seg.resize(4*n);
        lazy.resize(4*n, 0);
        arr = a; // if TLE use global array;
    }

    void build(int root, int low, int high){
        if(low > high) return;
        if(low == high){
            seg[root] = arr[low];
            return;
        }
        int mid = low + (high-low)/2;
        build(2*root+1, low, mid);
        build(2*root+2, mid+1, high);
        seg[root] = seg[root*2+1] + seg[2*root+2];
    }

    void pointUpdate(int root, int low, int high, int idx, int val){
        if(low == high){
            seg[idx] = val;
            return;
        }

        int mid = low + (high-low)/2;

        if(idx<=mid){
            pointUpdate(2*root+1, low, mid, idx, val);
        }else{
            pointUpdate(2*root+2, mid+1, high, idx, val);
        }

        seg[root] = seg[2*root+1] + seg[2*root+2];
    }

    int rangeQuery(int root, int low, int high, int l, int r){
        if( r<low ||  high < l || l > r)return 0;
        // full overlap
        if(l <= low && high <= r){
            return seg[root];
        }
        //partial overlap
        int mid = low + (high-low)/2;
        return rangeQuery(2*root+1, low, mid, l, r) + 
                rangeQuery(2*root+2, mid+1, high, l, r);
    }

    int rangeUpdate(int root, int low, int high, int l, int r, int val){
        if(lazy[root]!=0){
            seg[root]+=(high-low+1)*lazy[root];

            if(low!=high){
                lazy[2*root+1] += val;
                lazy[2*root+2] += val;
            }
            lazy[root] = 0;
        }

        if(l > high || r < low) return;

        if(low>=l && high<=r){
            seg[root] += (high-low+1) * val;
            if(low!=high){
                lazy[2*root+1]+=val;
                lazy[2*root+2]+=val;
            }
            return;
        }
        int mid  = low+(high-low)/2;
        rangeUpdate(2*root+1, low, mid, l, r, val);
        rangeUpdate(2*root+2, mid+1, high, l, r, val);
    }
};
```
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
## Binary lifting
```
vector<vector<int>> table;
int mx = 19; // this mx is quite large , any will be reahed only in few cases 
// so improve time we can set mx baaed on the value of N.
// so this max = ceil(LOG2(N)); 

void build(int n){
    table[0] = parent;

    for(int i = mx-1;i>=0;i--){
        for(int node = 1;node<=n;node++){
            table[i][node]= table[i-1][table[i-1][node]];
        }
    }
}

int calJump(int c){
    int curr = jump[c];
    if(occ[curr])return 0;
    int jmp = 1;
    for(int j = mx-1;j>=0;j--){
        int jp = table[j][curr]; //2^j th parent of curr;
        if(occ[jp]) continue;
        else{
            curr = jp;
            jmp += 1<<j;
        }
    }
    occ[curr] = true;
    return jmp;
}
void vivek()
{   
    int n;
    cin >> n;
    // one based indexing
    parent.resize(n+1);
    occ.resize(n+1, false);
    occ[0]= true;
    parent[1] = 0;
    parent[0] = 0;
    loop2(i, n){
        cin >> parent[i];
    }
     jump.resize(n+1);
     loop2(i, n){
        cin >> jump[i];
     } 

     table.resize(mx, vector<int>(n+1, 1));

     build(n);

    for(int i=1;i<=n;i++){
        int ans = calJump(i);
        cout<<ans<<endl;
    }
}

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
	long long res = 1;
	long long base = x % m;

	while (n > 0) {
		if (n % 2 == 1) {
			res = (res * base) % m;
		}
		base = (base * base) % m;
		n /= 2;
	}
	return (int)res;
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
## Find Factorial (genrally precomputed when required)
```

```
## Binary Exponentiation
```
long long power(long long a, long long b) {
    long long result = 1;
    while(b) {
        if (b & 1) 
        result = result * a;
        a = a * a;
        b >>= 1;
    }
    return result; //TIME : log (b);
}
```
## apply MOD when there can be overflow
```
ans = (ans + (1LL * arr[i] * (i - left) * (right - i)) % mod) % mod;
```

## Modular nCr
```
int modularnCr(int n , int r, int MOD){
	if(r<0 || r>n){
		return 0;

	}
	long long a = fact(n) % MOD;
	long long b = (fact(r)*fact(n-r)) % MOD;

	return a * power(b, MOD-2) % MOD;
}
```

## Mnacher Algo - Longest Pallingdromic Substring in O(n)
```
string findLongestPalindromicString(string text) {
        // using manacher Algo
        int N = text.length();
        if (N == 0)
            return "";
        if (N == 1)
            return text;

        N = 2 * N + 1;

        int L[N];
        L[0] = 0;
        L[1] = 1;

        int C = 1;

        int R = 2;

        int i = 0;

        int iMirror;
        int expand = -1;
        int diff = -1;
        int maxLPSLength = 0;
        int maxLPSCenterPosition = 0;
        int start = -1;
        int end = -1;

        for (i = 2; i < N; i++) {
            iMirror = 2 * C - i;
            expand = 0;
            diff = R - i;
            if (diff >= 0) {
                if (L[iMirror] < diff)
                    L[i] = L[iMirror];
                else if (L[iMirror] == diff && R == N - 1)
                    L[i] = L[iMirror];
                else if (L[iMirror] == diff && R < N - 1) {
                    L[i] = L[iMirror];
                    expand = 1;
                } else if (L[iMirror] > diff) {
                    L[i] = diff;
                    expand = 1;
                }
            } else {
                L[i] = 0;
                expand = 1;
            }

            if (expand == 1) {
                while (
                    ((i + L[i]) < N && (i - L[i]) > 0) &&
                    (((i + L[i] + 1) % 2 == 0) ||
                     (text[(i + L[i] + 1) / 2] == text[(i - L[i] - 1) / 2]))) {
                    L[i]++;
                }
            }
            if (L[i] > maxLPSLength) {
                maxLPSLength = L[i];
                maxLPSCenterPosition = i;
            }
            if (i + L[i] > R) {
                C = i;
                R = i + L[i];
            }
        }

        start = (maxLPSCenterPosition - maxLPSLength) / 2;
        end = start + maxLPSLength - 1;
        // cout << "LPS of string is " << text << " : ";
        string ans;
        for (i = start; i <= end; i++)
            ans.push_back(text[i]);

        return ans;
    }
```

## sample DIGIT DP code 
```
#include <bits/stdc++.h>
#define int long long int
#define no cout<<"NO"<<endl;
#define yes cout<<"YES"<<endl;
#define Sort(v) sort(v.begin(), v.end())
#define Sortr(v) sort(v.rbegin(), v.rend())
#define vecpair vector<pair<int, int>> 
#define loop1(i,n) for(int i=0;i<n;i++)
#define loop2(i,n) for(int i=1;i<=n;i++)
#define print1(x) cout<<x<<endl;
#define print2(x) cout<<x<<" ";
#define NL cout<<"\n";
#define umap unordered_map
#define all(v) v.begin(), v.end()
using namespace std;
// to print cout of digit 3 in all the numbers form l to r
// digit DP sort of template 
int dp[20][2][20]; //[index][tight][count]
int solve(string &s, int idx, int tight, int count){
    if(idx == s.length()) return count;
    if(dp[idx][tight][count]!=-1)return dp[idx][tight][count];
    int limit = (tight == 1 ? s[idx]-'0' : 9);
    int ans = 0;
    for(int i=0;i<=limit;i++){
        int update = count+(i==3)?1:0;
        ans+=solve(s, idx+1, (tight & (i==s[idx]-'0')), update);
    }
    return dp[idx][tight][count] = ans;
}
void vivek()
{   
    int l, r;
    cin >> l >> r;
    
    string ri = to_string(r);
    memset(dp, -1, sizeof(dp));
    int rightAns = solve(ri, 0, 1, 0); //idx, tight=1, count = 0 -> tight = 1 initiall since once false it remail false;
    string le = to_string(l-1);
    memset(dp, -1, sizeof(dp));
    int leftAns = solve(le, 0, 1, 0);
    
    cout<< rightAns - leftAns<<endl;
}

signed main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int t = 1;
    // /*is Single Test case?*/ cin >> t;
    while (t--)
    {
        vivek();
    }

    cerr << "time taken : " << (float)clock() / CLOCKS_PER_SEC << " secs" << endl;
    return 0;
}

```
