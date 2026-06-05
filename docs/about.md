# About this project

The original purpose of the project was exploratory: walk through the steps in the DADA2 workflow to understand the underlying implementation for each key command, initially to replicate results from the R DADA2 workflow but to also explore potential paths for improving the implementation, as we use this quite extensively in our own research. 

The current implementation:
- fully implements the DADA2 workflow, 
- closely matches results from DADA2, and 
- runs faster than the original implementatin for Illumina and PacBio data, using less memory

See the latest [benchmarking results](benchmarking.md) for the most up-to-date information.

---

## Learning

This was also meant to be a learning opportunity on several fronts. I have been learning Rust over the last few years in my (vanishingly small) spare time, and actually planned a general port of DADA2 a few years back. I had programmed various languages over the years (Perl, Python, R, C++, and a little Java), but apart from R and Python I'm pretty, um, rusty.

## AI

I have noticed a dramatic improvement in coding-based AI tools and agents, also noted by many others in the field. However there has been a lot of discussion on the best approaches to reimplementation strategies involving these approaches, and many controversies along the way. 

I've long been involved in open-source development and open-science projects, and I'm also the leader of a bioinformatics core group. I understand the community responses to this as well as the potential benefits to this community, and I believe community standards are needed that ensure we follow some general guidelines that ensure some degree of consistency and support, provenance that ensure the original implementors retain clear credit for their work, and where community members can play an active role. 

The reality is: many open-source projects, including a few I am directly involved in, are struggling to reconcile with the changes occurring and how best to responsibly utilize these tools, as well as how to accept contributions that also use them.

## Guidelines

AI is having a clear, fundamental, and disruptive impact, and simply ignoring this is to one's detriment. Some of us are also mentors, and as such we have an obligation to understand this wildly changing landscape and help prepare students, postdocs, and scientist on how to best use AI for their future career path. However many controveries exist over the use of these tools, in some cases arguably crossing moral and ethical boundaries. Standards are sorely needed.

Thankfully, within the bioinformatics community these are starting to coalesce, for example [rewrites.bio](https://rewrites.bio). 

Therefore, this project will follow [rewrites.bio](https://rewrites.bio) guidelines as closely as possible to: 

* Ensure the original DADA2 developers and contributors are acknowledged,
* Follows the original implementation's details,
* Include tests and benchmarks for the work,
* Utilize consistent libraries where possible,
* Release as open source and follow the original licensing
* Should interest arise: develop a community that can contribute.

## Caveat

One clear caveat where this implementation will vary from the rewrites.bio standards: due to key implementation details (conversion of R/C++ to Rust including error models, a fully ported R-based LOESS implementation in Rust), results will vary slightly. However we strive to reproduce results as closely as possible, within reason, and outline where we see differences so you can make your own judgement. We have also added the ability to use R and Python for custom error model analysis to more closely emulate what we see from the original implementation.

We also do not want to prevent additional outcomes or functionality that may come from the work in this project by being constrained to solely emulating the original code. For example a key outcome has come from exposing the k-mer size as an option: **alternative, longer k-mer lengths improve performance for PacBio denoising**. This is something that needs to be explored more closely, but appears to result in a substantial improvement in data processing with no apparent difference in ASV results.