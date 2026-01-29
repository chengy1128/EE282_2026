# Homework 2: Practice Problems and Answer Key

## Problem 1: Navigating and manipulating directories (`cd`, `../`, `mkdir`, `rmdir`)

You are currently located in the directory:

- `~/teaching/labs/week2`

Assume the following is true:

- `~/teaching/labs/week2` already exists.
- `~/teaching/labs/week3` does **not** exist yet.
- There is **no** directory named `scratch` anywhere in `~/teaching/labs/`.

Write a sequence of shell commands that accomplishes the tasks below **in order**:

1. Move to the parent directory of `week2` (i.e., into `~/teaching/labs/`) using `cd` and `../`.
2. Create a new directory named `week3` inside `~/teaching/labs/`.
3. Move into `week3`.
4. Create a subdirectory inside `week3` named `scratch`.
5. Move into `scratch`.
6. Move back up to `week3` using `cd` and `../`.
7. Remove the (now-empty) `scratch` directory using `rmdir`.
8. End in the `~/teaching/labs/` directory.

### Answer key (Problem 1)

```bash
# 1) go to the parent of week2 (i.e., labs)
cd ../

# 2) create week3 inside labs
mkdir week3

# 3) move into week3
cd week3

# 4) create scratch inside week3
mkdir scratch

# 5) move into scratch
cd scratch

# 6) back up to week3
cd ../

# 7) remove the empty scratch directory
rmdir scratch

# 8) end in ~/teaching/labs/
cd ../
```


## Problem 2: R indexing with text column names (matrix vs. data frame)

Consider the following R code:

```r
mymatrix <- matrix(
  1:6,
  nrow = 3,
  ncol = 2,
  byrow = FALSE,
  dimnames = list(NULL, c("col1", "col2"))
)

mydf <- data.frame(
  col1 = 1:3,
  col2 = 4:6
)
```

For **each** of the following expressions, state:

1. The values it returns (what you would see printed in the Console).
2. The **class** of the returned object.
3. Its **shape**:
   - for vectors: `length`
   - for matrices/data frames: `dim`

Expressions to evaluate:

- `mymatrix[,'col1']`
- `mydf[,'col1']`
- `mydf['col1']`
- `mydf$col1`
- `mydf[['col1']]`

### Answer key (Problem 2)

The key idea:

- A **matrix** is a 2D atomic object; `mymatrix[ , 'col1']` selects a column and (by default) **drops** to a vector.
- A **data frame** is a list of columns; different indexing operators (`[`, `[[`, `$`) return different types.

| Expression | Printed value | Class | Shape |
|---|---|---|---|
| `mymatrix[,'col1']` | `1 2 3` | `integer` | `length = 3` |
| `mydf[,'col1']` | `1 2 3` | `integer` | `length = 3` |
| `mydf['col1']` | a 1-column data frame with `col1` | `data.frame` | `dim = c(3, 1)` |
| `mydf$col1` | `1 2 3` | `integer` | `length = 3` |
| `mydf[['col1']]` | `1 2 3` | `integer` | `length = 3` |

More explicit printed output examples:

```r
mymatrix[,'col1']
# [1] 1 2 3

mydf[,'col1']
# [1] 1 2 3

mydf['col1']
#   col1
# 1    1
# 2    2
# 3    3

mydf$col1
# [1] 1 2 3

mydf[['col1']]
# [1] 1 2 3
```

Notes on *why* they differ:

- `mydf['col1']` uses single-bracket indexing on a data frame, which returns a **data frame** (a subset of columns).
- `mydf[['col1']]` and `mydf$col1` extract the **column vector itself**.
- `mydf[,'col1']` looks like matrix-style 2D indexing; by default it also returns the extracted column vector (it *drops* to a vector). If you wanted a 1-column data frame instead, you would use:
  - `mydf[, 'col1', drop = FALSE]`


## Problem 3: Sharing access using symbolic permissions (`chmod u/g/o +/- r/w/x`)

You are user `alex`. Your colleague is **not** in your Unix group, so you must treat them as “other” (`o`) for permissions.

You have this directory and file inside your home directory:

- Directory: `~/collab`
- File: `~/collab/notes.txt`

Current permissions are:

```bash
$ ls -ld ~ ~/collab ~/collab/notes.txt
drwx------  1 alex alex 4096 Jan 29 12:00 /home/alex
drwx------  1 alex alex 4096 Jan 29 12:01 /home/alex/collab
-rw-------  1 alex alex  210 Jan 29 12:02 /home/alex/collab/notes.txt
```

Goal:

- Your colleague should be able to:
  1. Traverse into `~/collab` (so `cd /home/alex/collab` works).
  2. List the contents of `~/collab` (so `ls /home/alex/collab` works).
  3. Read the file (so `cat /home/alex/collab/notes.txt` works).
- Your colleague should **not** be able to list your home directory (so `ls /home/alex` should still fail due to lack of read permission).

Write the `chmod` command(s) (symbolic form only) needed to meet these requirements.

### Answer key (Problem 3)

```bash
# Allow others to *traverse* your home directory, but not list it.
chmod o+x ~

# Allow others to traverse and list the collab directory.
chmod o+rx ~/collab

# Allow others to read the file.
chmod o+r ~/collab/notes.txt
```

What the permissions become:

```bash
$ ls -ld ~ ~/collab ~/collab/notes.txt
drwx-----x  1 alex alex 4096 Jan 29 12:00 /home/alex
drwx---r-x  1 alex alex 4096 Jan 29 12:01 /home/alex/collab
-rw----r--  1 alex alex  210 Jan 29 12:02 /home/alex/collab/notes.txt
```

Why this works:

- The colleague needs **execute** (`x`) on every parent directory to reach the file (here: `~` and `~/collab`).
- The colleague needs **read** (`r`) on `~/collab` to run `ls` there.
- The colleague needs **read** (`r`) on the file itself to view its contents.
- Keeping `~` without `o+r` prevents listing the entire home directory.
