# sb_4.11

Rename: ```SeqBiome__STC``` 

Or possibly ```SB__STC```

Even more likely, fork this and make it open.


quickref:
```
# SSH clone 
git clone git@github.com:handibles/sb_4.11
# HTTPS clone
git clone https://github.com/handibles/sb_4.11


# commit
git init
git add .
git commit -m "first stash and commit"                                 # note comment arg
git remote add origin github@github.com:handibles/XXX_REPO_NAME_XXX    # note syntax = ssh
git push -u origin master                                              # local -> cloud
```

This repo should be pullable and usable for SB projects (16S focused at the mo).

general structure

```
      analysis
        ¬| bash - pipeline, rsync, misc.
         | py   - misc.
         | R    - scripts and functions

      documents
        ¬| reproducibility guidelines etc
         | outputs, generated HTML or PDF docs

      input
        ¬| external data

      output
        ¬| outputs from different stages 

      vis
        ¬| visual outputs

      SB__4.11.Rproj    # note prob with nomenclature

```

See 0.0 for the nascent blog-ish post detailing creating this repo. Should ideally be moved to a more narrative repo of unknown descent. 



