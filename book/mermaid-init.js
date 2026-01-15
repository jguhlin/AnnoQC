(function () {
  function renderMermaidBlocks() {
    if (typeof mermaid === "undefined") return;

    // Convert fenced code blocks (```mermaid) into <div class="mermaid"> nodes.
    const blocks = document.querySelectorAll(
      "pre > code.language-mermaid, pre > code.mermaid"
    );
    for (const code of blocks) {
      const pre = code.parentElement;
      if (!pre) continue;

      const wrapper = document.createElement("div");
      wrapper.className = "mermaid-wrapper";
      const container = document.createElement("div");
      container.className = "mermaid";
      container.textContent = code.textContent || "";
      wrapper.appendChild(container);
      pre.replaceWith(wrapper);
    }

    // Mermaid v9 API
    mermaid.initialize({
      startOnLoad: false,
      flowchart: { useMaxWidth: false },
      themeVariables: { fontSize: "11px" },
    });
    mermaid.init(undefined, document.querySelectorAll(".mermaid"));

    // Post-process: Mermaid often emits <svg width="100%"> which prevents horizontal
    // scrolling and can make wide diagrams appear clipped. Force an explicit pixel
    // width based on the viewBox so wrappers can scroll.
    for (const svg of document.querySelectorAll(".mermaid > svg")) {
      const viewBox = svg.getAttribute("viewBox");
      if (!viewBox) continue;
      const parts = viewBox.trim().split(/\s+/);
      if (parts.length !== 4) continue;
      const vbWidth = Number(parts[2]);
      if (!Number.isFinite(vbWidth) || vbWidth <= 0) continue;
      svg.style.width = `${Math.ceil(vbWidth)}px`;
      svg.style.maxWidth = "none";
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", renderMermaidBlocks);
  } else {
    renderMermaidBlocks();
  }
})();
