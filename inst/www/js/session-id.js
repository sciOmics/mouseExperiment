// session-id.js — CODE_REVIEW.md F.7, layer 3 (cookie setter)
//
// On first connect, check for the mE_session cookie. If absent,
// generate a 256-bit random session ID via crypto.getRandomValues and
// set the cookie. Subsequent connects (including reconnects after
// disconnect) reuse the existing cookie.
//
// The server NEVER reads the cookie from JS — it reads the HTTP Cookie
// header directly on each request. This file only exists to *set* the
// cookie the first time a user lands on the dashboard. Once set, the
// browser handles sending it on every subsequent HTTP request to the
// same origin.
//
// Cookie attributes:
//   SameSite=Strict — no cross-site request can submit the cookie
//   Secure          — HTTPS only (browser refuses to set on http://)
//   Path=/          — scoped to this origin's app
//   max-age=86400   — 24-hour expiry; user re-identifies after a day
//
// Not HttpOnly — JS sets it client-side, so the cookie is readable
// from JS. XSS risk is mitigated by auditing user-supplied strings in
// the dashboard's HTML output.

(function() {
  function hasSessionCookie() {
    return document.cookie.split(";").some(function(c) {
      return c.trim().indexOf("mE_session=") === 0;
    });
  }

  function generateSid() {
    // 32 bytes = 256 bits of CSPRNG output, hex-encoded to 64 chars.
    var buf = new Uint8Array(32);
    crypto.getRandomValues(buf);
    return Array.from(buf, function(b) {
      return ("0" + b.toString(16)).slice(-2);
    }).join("");
  }

  function setSessionCookie() {
    if (hasSessionCookie()) return;
    var sid = generateSid();
    // Build cookie attribute string. Use Secure only on HTTPS; on
    // localhost dev (http://) the browser would refuse Secure cookies.
    var attrs = [
      "mE_session=" + sid,
      "Path=/",
      "Max-Age=86400",
      "SameSite=Strict"
    ];
    if (window.location.protocol === "https:") {
      attrs.push("Secure");
    }
    document.cookie = attrs.join("; ");
  }

  // Run synchronously so the cookie exists before Shiny opens its
  // websocket. If we deferred to DOMContentLoaded, the first server
  // tick could see a cookieless request and the banner observer would
  // briefly miss the user's jobs on reconnect.
  setSessionCookie();
})();
