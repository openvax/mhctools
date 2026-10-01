"""Check relative links and anchors in the documentation.

Run ``python scripts/check_docs_links.py``; it exits 1 and lists every broken
link. ``tests/test_docs_links.py`` runs the same check. Links from ``docs/``
must stay inside ``docs/`` (use a GitHub URL for repository files), because the
site is built from that directory only.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
SITE = "https://openvax.github.io/mhctools/"

LINK = re.compile(r"(?<!\!)\[(?:[^\]]|\n)*?\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$", re.M)
FENCE = re.compile(r"^```.*?^```", re.M | re.S)


def slug(heading):
    """GitHub/MkDocs anchor for a heading."""
    text = re.sub(r"`([^`]*)`", r"\1", heading)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = text.lower().strip()
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def anchors(path):
    text = FENCE.sub("", path.read_text())
    found = {}
    for _, heading in HEADING.findall(text):
        base = slug(heading)
        n = found.get(base, -1) + 1
        found[base] = n
        if n:
            found["%s-%d" % (base, n)] = 0
    return set(found)


def _resolve_site_url(url):
    """Map a site URL to the docs file it will serve, or None."""
    rest = url[len(SITE):].strip("/")
    if not rest:
        return DOCS / "index.md"
    for candidate in (DOCS / (rest + ".md"), DOCS / rest / "index.md"):
        if candidate.exists():
            return candidate
    return None


def check():
    errors = []
    files = sorted(DOCS.rglob("*.md")) + [ROOT / "README.md"]
    cache = {}

    def anchors_of(path):
        if path not in cache:
            cache[path] = anchors(path)
        return cache[path]

    for source in files:
        text = FENCE.sub("", source.read_text())
        rel_source = source.relative_to(ROOT)
        in_docs = DOCS in source.parents
        for target in LINK.findall(text):
            if target.startswith("mailto:"):
                continue
            if target.startswith("#"):
                if target[1:] not in anchors_of(source):
                    errors.append("%s: missing anchor %s" % (rel_source, target))
                continue
            if target.startswith(("http://", "https://")):
                if not target.startswith(SITE):
                    continue
                path_part, _, anchor = target.partition("#")
                resolved = _resolve_site_url(path_part)
                if resolved is None:
                    errors.append("%s: site link has no page: %s" % (rel_source, target))
                elif anchor and anchor not in anchors_of(resolved):
                    errors.append("%s: missing anchor in %s: #%s" % (
                        rel_source, resolved.relative_to(ROOT), anchor))
                continue
            path_part, _, anchor = target.partition("#")
            resolved = (source.parent / path_part).resolve()
            if in_docs and DOCS not in resolved.parents and resolved != DOCS:
                errors.append(
                    "%s: link leaves docs/ (use a GitHub URL): %s" % (rel_source, target))
            elif not resolved.exists():
                errors.append("%s: missing file: %s" % (rel_source, target))
            elif anchor and resolved.suffix == ".md" and anchor not in anchors_of(resolved):
                errors.append("%s: missing anchor in %s: #%s" % (
                    rel_source, resolved.relative_to(ROOT), anchor))
    return errors


def main():
    errors = check()
    for error in errors:
        print(error)
    print("%d broken link(s)" % len(errors))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
